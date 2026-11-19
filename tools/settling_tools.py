"""
settling_tools: settling-velocity fits for a stirred mixing tank from an OBS decay record
(single size or size distribution), turbulence calibration, and calibration files.
(Renamed from settling_single.py in the calibration-file release.)

Model (bounded tank, 1-D vertical, z measured from the surface in the solver):
    dc/dt + d/dz( w c - K(t) dc/dz ) = 0,    K(t) = Kmin + K0 exp(-t/tau)
    * zero particle flux at the free surface
    * depositional bed (particles that reach the bed stay there)
    * initial state = steady stirred profile for K0 (nearly uniform when w h / K0 is small)
Model output is normalised by its own initial reading at the sensor, exactly as the data are
normalised by the mixer-on plateau:  c = (NTU - Cbackground) / (Co - Cbackground).

Inputs you control:  particle size, geometry, water T/S, K0 / tau (from the ADV), settling-law coefficients.
Fitted:              the settling velocity w (optionally K0 and tau too, as a diagnostic).
Afterwards:          invert_ws() turns w into particle density (nf fixed) or fractal dimension (density fixed).

Turbulence can be given either as K(t) = Kmin + K0 exp(-t/tau) or, new in v03, as an ADV-derived shape:
K(t) = Kmin + K0 * shape(t), with shape(t) from adv_shape() (shape(0) = 1, so K0 is the stirred-state diffusivity).

Geometry convention follows SettlingFit_OBS_PolySize: water depth h and sensor height above the bottom z_obs.
"""
import json
import datetime
from dataclasses import dataclass
import numpy as np
import pandas as pd
import warnings
from scipy.optimize import least_squares, brentq
from settling_solver import solve_single, solve_classes, K_exp

# finite-difference step (relative to the log-parameters) for the least-squares Jacobians; the ODE solver output carries
# ~1e-5 integration noise, so the scipy default (~1e-8) gives unreliable derivatives and uncertainties
DIFF_STEP = 1e-3


# ----------------------------------------------------------------------------- water and settling law
def density(To, Sppt):
    rho_fresh = 1000 * (1 - (To + 288.9414) / (508929.2 * (To + 68.12963)) * (To - 3.9863) ** 2)
    Acoef = 0.824493 - 0.0040899 * To + 0.000076438 * To ** 2 - 0.00000082467 * To ** 3 + 0.0000000053675 * To ** 4
    Bcoef = -0.005724 + 0.00010227 * To - 0.0000016546 * To ** 2
    return rho_fresh + Acoef * Sppt + Bcoef * Sppt ** (3 / 2) + 0.00048314 * Sppt ** 2


def visc(To):
    return 2.7488e-07 + 1.4907e-06 * np.exp(-0.034812 * To)          # kinematic viscosity [m2/s]


COEFS = {1: (20, 1.258),   # flocs matching Stokes at the primary particles
         2: (120, 0.0),    # viscous flocs
         3: (18, 0.866),   # sieve-sized sand / smooth spheres (Ferguson & Church C1=18, C2=1)
         4: (20, 0.91),    # nominal diameter natural sand
         5: (18, 0.548)}   # sphere


def settlingvel(d_um, dp_um, nf, coef, rhos, T, S):
    """Settling velocity [m/s] (Ferguson & Church adapted to fractal flocs); nf=3 gives a solid particle."""
    rho, nu, g = density(T, S), visc(T), 9.81
    Rs = (rhos - rho) / rho
    b1, b2 = COEFS.get(coef, COEFS[5])
    d, dp = np.asarray(d_um, float) / 1e6, dp_um / 1e6
    return (g * Rs * d ** (nf - 1)) / (b1 * nu * dp ** (nf - 3) + b2 * np.sqrt(g * Rs * d ** nf * dp ** (nf - 3)))


def invert_ws(ws, d_um, which, T, S, coef, rhos=None, nf=None, dp_um=8.0):
    """Back out particle density (which='rhos', needs nf; solid particle: nf=3) or
    fractal dimension (which='nf', needs rhos and dp_um) from a settling velocity. NaN if no solution."""
    rho_w = density(T, S)
    try:
        if which == "rhos":
            f = lambda r: settlingvel(d_um, dp_um, nf, coef, r, T, S) - ws
            return brentq(f, rho_w + 1e-3, 6000.0)
        if which == "nf":
            f = lambda n: settlingvel(d_um, dp_um, n, coef, rhos, T, S) - ws
            return brentq(f, 1.2, 3.0)
    except ValueError:
        return np.nan
    raise ValueError("which must be 'rhos' or 'nf'")


# ----------------------------------------------------------------------------- data
def load_decay(path, t_col="sec", ntu_col="NTU", co_col="Co_NTU_avg", cb_col="Cbackground_NTU_avg"):
    """Read the *_decayONLY.csv. Co and Cbackground are the first-row values. Returns a dict with
    t [s since mixer off], ntu, Co, Cb and the normalised concentration c = (NTU-Cb)/(Co-Cb)."""
    df = pd.read_csv(path)
    t, ntu = df[t_col].to_numpy(float), df[ntu_col].to_numpy(float)
    Co, Cb = float(df[co_col].iloc[0]), float(df[cb_col].iloc[0])
    return dict(t=t, ntu=ntu, Co=Co, Cb=Cb, c=(ntu - Cb) / (Co - Cb))


def t50_estimate(t, c, window=7):
    """First time the (running-mean smoothed) normalised concentration drops below 0.5, by interpolation."""
    cs = pd.Series(c).rolling(window, center=True, min_periods=1).mean().to_numpy()
    i = int(np.argmax(cs <= 0.5))
    if i == 0:
        return np.nan
    return float(np.interp(0.5, [cs[i], cs[i - 1]], [t[i], t[i - 1]]))


@dataclass
class Geometry:
    h: float        # water depth [m]
    z_obs: float    # sensor height above the bottom [m]
    @property
    def z_s(self):  # sensor depth below the free surface [m] (what the solver uses)
        return self.h - self.z_obs
    @property
    def fall(self): # distance a particle must fall to clear the sensor [m]
        return self.h - self.z_obs


# ----------------------------------------------------------------------------- forward model and fit
def model_c(t, ws, K0, tau, Kmin, geom, N=600, Kfun=None):
    """Normalised sensor concentration at times t [s since mixer off] for a single class.
    Kfun (callable K(t)) overrides K0, tau, Kmin when given (e.g. from shape_K)."""
    t = np.asarray(t, float)
    has0 = t[0] <= 0.0
    teval = t if has0 else np.r_[0.0, t]
    with warnings.catch_warnings():
        # scipy BDF's sparse finite-difference Jacobian can overflow its step factor on flat regions; harmless here
        warnings.filterwarnings('ignore', category=RuntimeWarning, module='scipy.integrate')
        cm = solve_single(ws, Kfun if Kfun is not None else K_exp(K0, tau, Kmin), geom.h, geom.z_s, teval,
                          N=N, rtol=1e-5, atol=1e-8)
    return cm / cm[0] if has0 else cm[1:] / cm[0]


def fit_ws(t, c, geom, K0, tau, Kmin=5e-8, ws0=None, t_start=0.0, t_end=None, fit_turbulence=False, N=400,
           shape=None):
    """Least-squares fit of the settling velocity (and optionally K0, tau). Returns a dict.
    If `shape` (from adv_shape) is given, K(t) = Kmin + K0*shape(t), tau is unused and fit_turbulence fits K0 only.
    ws_sigma is the Jacobian-based 1-sigma, inflated for lag-1 autocorrelation of the residuals."""
    t, c = np.asarray(t, float), np.asarray(c, float)
    m = (t >= t_start) & ((t <= t_end) if t_end is not None else np.ones_like(t, bool))
    tt, cc = t[m], c[m]
    ws0 = ws0 if ws0 is not None else 7e-4
    nturb = 0 if not fit_turbulence else (1 if shape is not None else 2)
    def unpack(x):
        if nturb == 0:
            return np.exp(x[0]), K0, tau
        return (np.exp(x[0]), np.exp(x[1]), np.exp(x[2]) if nturb == 2 else tau)
    def resid(x):
        w, K, ta = unpack(x)
        Kf = shape_K(shape, K, Kmin) if shape is not None else None
        return model_c(tt, w, K, ta, Kmin, geom, N, Kfun=Kf) - cc
    x0 = [np.log(ws0)] + ([np.log(K0), np.log(tau)][:nturb])
    lo = [np.log(1e-5)] + ([np.log(1e-7), np.log(1.0)][:nturb])
    hi = [np.log(5e-2)] + ([np.log(1e-1), np.log(900.0)][:nturb])
    r = least_squares(resid, x0, bounds=(lo, hi), diff_step=DIFF_STEP)
    n, p = len(cc), len(x0)
    rho1 = np.corrcoef(r.fun[:-1], r.fun[1:])[0, 1]
    infl = (1 + max(rho1, 0)) / (1 - min(max(rho1, 0), 0.95))
    cov = np.linalg.pinv(r.jac.T @ r.jac) * np.sum(r.fun ** 2) / max(n - p, 1) * infl
    w, K, ta = unpack(r.x)
    out = dict(ws=w, K0=K, tau=ta, ws_sigma=w * np.sqrt(cov[0, 0]), rms=float(np.sqrt(np.mean(r.fun ** 2))),
               rho1=float(rho1), n=n, t=tt, c=cc, resid=r.fun, success=r.success)
    return out


def sensitivity_table(t, c, geom, K0_list, tau_list, d_um, T, S, coef, nf=3.0, dp_um=8.0, Kmin=5e-8,
                      ws0=None, t_start=0.0, t_end=None, N=200, shape=None):
    """Refit w for every (K0, tau) pair; shows how much the turbulence assumption moves the answer.
    With `shape` given the ADV shape is used and only K0_list matters (tau column is NaN)."""
    rows = []
    for K0 in K0_list:
        for tau in ([np.nan] if shape is not None else tau_list):
            f = fit_ws(t, c, geom, K0, 30.0 if shape is not None else tau, Kmin, ws0, t_start, t_end, N=N, shape=shape)
            rows.append(dict(K0=K0, tau=tau, ws_mm_s=f["ws"] * 1e3, rms=f["rms"],
                             rhos=invert_ws(f["ws"], d_um, "rhos", T, S, coef, nf=nf, dp_um=dp_um)))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------- calibrate turbulence with a standard
def calibrate_turbulence(t, c, geom, ws_known, free=("K0", "tau"), K0=1e-4, tau=10.0, Kmin=5e-8,
                         t_start=0.0, t_end=None, N=300, bounds=None):
    """Use a standard of known settling velocity (known size and density) to calibrate the turbulence:
    w is FIXED at ws_known and the parameters named in `free` (any of 'K0', 'tau', 'Kmin') are fitted to the
    decay record.  Returns fitted values, rms, a lag-1-autocorrelation-inflated 1-sigma for each free
    parameter and the residuals.  Parameters not in `free` stay at the values passed in."""
    t, c = np.asarray(t, float), np.asarray(c, float)
    m = (t >= t_start) & ((t <= t_end) if t_end is not None else np.ones_like(t, bool))
    tt, cc = t[m], c[m]
    start = dict(K0=K0, tau=tau, Kmin=Kmin)
    bnd = dict(K0=(1e-7, 1e-1), tau=(0.5, 900.0), Kmin=(5e-9, 1e-2))
    if bounds:
        bnd.update(bounds)
    free = list(free)
    def unpack(x):
        p = dict(start)
        p.update({k: float(np.exp(v)) for k, v in zip(free, x)})
        return p
    def resid(x):
        p = unpack(x)
        return model_c(tt, ws_known, p["K0"], p["tau"], p["Kmin"], geom, N) - cc
    r = least_squares(resid, [np.log(start[k]) for k in free],
                      bounds=([np.log(bnd[k][0]) for k in free], [np.log(bnd[k][1]) for k in free]), diff_step=DIFF_STEP)
    rho1 = np.corrcoef(r.fun[:-1], r.fun[1:])[0, 1]
    infl = (1 + max(rho1, 0)) / (1 - min(max(rho1, 0), 0.95))
    cov = np.linalg.pinv(r.jac.T @ r.jac) * np.sum(r.fun ** 2) / max(len(cc) - len(free), 1) * infl
    p = unpack(r.x)
    out = dict(p, rms=float(np.sqrt(np.mean(r.fun ** 2))), rho1=float(rho1), free=free, resid=r.fun, t=tt, c=cc)
    out["sigma_ln"] = dict(zip(free, np.sqrt(np.diag(cov))))      # 1-sigma in ln(parameter), i.e. fractional
    return out



# ----------------------------------------------------------------------------- ADV-derived turbulence shape
def adv_shape(path, t_col="decay_time_s", vel_cols=("v1_m_s", "v2_m_s", "v3_m_s"), window_s=2.0,
              noise=None, noise_from_s=60.0, t_cut_s=30.0, smooth=3, normalize="first"):
    """Dimensionless shape of the turbulence decay from an ADV record (K(t) = K0 * shape(t) + Kmin).

    Steps: (1) windows of `window_s` seconds (50 % overlap); in each the linear trend of every velocity
    component is removed and the fluctuation variance v1'^2+v2'^2+v3'^2 is averaged (this strips the slowly
    varying mean flow and offsets, which the raw sum of squares contains); (2) the noise floor (median of that
    quantity for t > `noise_from_s`, or the value passed as `noise`) is subtracted; (3) K is taken proportional
    to the rms fluctuation velocity, i.e. shape = sqrt(excess); (4) a `smooth`-window running mean is applied,
    the shape is set to zero for t > `t_cut_s` (indistinguishable from noise; the persistent weak mixing is
    carried by Kmin instead) and normalised so shape(first window) = 1 (normalize="first") or so that its
    maximum is 1 (normalize="max", better if the record starts with a weak spin-up window).
    `noise` must be in the same units as the detrended fluctuation variance (m2/s2), which is not exactly the
    'background energy' of the raw sum of squares.
    Returns dict(t, shape, noise, fluct (variance per window), umax)."""
    d = pd.read_csv(path)
    ta = d[t_col].to_numpy(float)
    V = d[list(vel_cols)].to_numpy(float)
    dt = float(np.median(np.diff(ta)))
    n = max(int(round(window_s / dt)), 8)
    ts, fl = [], []
    for i in range(0, len(ta) - n, max(n // 2, 1)):
        seg = V[i:i + n]
        x = np.arange(n) - 0.5 * (n - 1)
        resid = seg - seg.mean(0) - np.outer(x, (x @ (seg - seg.mean(0))) / (x @ x))   # linear detrend
        ts.append(ta[i:i + n].mean()); fl.append((resid ** 2).mean(0).sum())
    ts, fl = np.array(ts), np.array(fl)
    if noise is None:
        late = ts > noise_from_s
        if late.sum() < 5:
            raise ValueError("ADV record too short to estimate the noise floor; pass noise= or lower noise_from_s")
        noise = float(np.median(fl[late]))
    u = pd.Series(np.sqrt(np.clip(fl - noise, 0.0, None))).rolling(smooth, center=True, min_periods=1).mean().to_numpy().copy()
    u[ts > t_cut_s] = 0.0
    ref = u.max() if normalize == "max" else u[0]
    if ref <= 0:
        raise ValueError("no turbulence above the noise floor in the first window")
    return dict(t=ts, shape=u / ref, noise=noise, fluct=fl, umax=float(ref))


def shape_K(shape, K0, Kmin=0.0):
    """Callable K(t) = K0 * shape(t) + Kmin; `shape` is the dict returned by adv_shape().
    Held at shape[0] before the first window centre, zero after the last."""
    ts, sh = shape["t"], shape["shape"]
    return lambda t: K0 * float(np.interp(t, ts, sh, left=sh[0], right=0.0)) + Kmin


# ----------------------------------------------------------------------------- size range (several classes)
def size_classes(d_lo, d_hi, n=7, weights=None):
    """Split a nominal size range into n equal-width bins; returns bin-centre diameters (um) and
    mass fractions (sum 1).  Default: uniform in diameter.  Pass `weights` (length n) to override,
    e.g. from a measured PSD."""
    edges = np.linspace(d_lo, d_hi, n + 1)
    d = 0.5 * (edges[1:] + edges[:-1])
    p = np.ones(n) / n if weights is None else np.asarray(weights, float) / np.sum(weights)
    return d, p


def class_ws(d_um, which_known, T, S, coef, rhos=None, nf=3.0, dp_um=8.0):
    """Settling velocity (m/s) of every class diameter for a known density (solid: nf=3)."""
    return np.array([settlingvel(d, dp_um, nf, coef, rhos, T, S) for d in np.atleast_1d(d_um)])


def model_c_classes(t, ws, p, K0, tau, Kmin, geom, N=300, Kfun=None):
    """Normalised sensor concentration for several classes (velocities ws, fractions p), one integration.
    Kfun (callable K(t)) overrides K0, tau, Kmin when given."""
    t = np.asarray(t, float)
    has0 = t[0] <= 0.0
    teval = t if has0 else np.r_[0.0, t]
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=RuntimeWarning, module='scipy.integrate')
        cm = solve_classes(ws, p, Kfun if Kfun is not None else K_exp(K0, tau, Kmin), geom.h, geom.z_s, teval, N=N)
    return cm / cm[0] if has0 else cm[1:] / cm[0]


def fit_ws_classes(t, c, geom, d_um, p, K0, tau, Kmin=5e-8, ws50_0=None, t_start=0.0, t_end=None, beta=2.0, N=250,
                   shape=None):
    """Fit one velocity scale for a size distribution: ws_i = ws50 * (d_i/d50)^beta (beta ~2 Stokes, ~1.3 for
    75-90 um sand-like spheres in this law; pass the log-slope of your settling law).  Independent of any
    density / fractal-dimension assumption.  Returns ws50, its 1-sigma (autocorrelation-inflated), the class
    velocities, rms, lag-1 residual autocorrelation and residuals."""
    t, c = np.asarray(t, float), np.asarray(c, float)
    m = (t >= t_start) & ((t <= t_end) if t_end is not None else np.ones_like(t, bool))
    d_um, p = np.asarray(d_um, float), np.asarray(p, float)
    d50 = _d50(d_um, p)
    dshape = (d_um / d50) ** beta
    Kf = shape_K(shape, K0, Kmin) if shape is not None else None
    ws50_0 = ws50_0 or 7e-4
    r = least_squares(lambda x: model_c_classes(t[m], np.exp(x[0]) * dshape, p, K0, tau, Kmin, geom, N, Kfun=Kf) - c[m],
                      [np.log(ws50_0)], bounds=([np.log(1e-5)], [np.log(5e-2)]), diff_step=DIFF_STEP)
    ws50 = float(np.exp(r.x[0]))
    rho1, sig = _sigma(r, 1)
    return dict(ws50=ws50, ws50_sigma=ws50 * sig[0], ws=ws50 * dshape, d50=d50, rms=float(np.sqrt(np.mean(r.fun ** 2))),
                rho1=rho1, resid=r.fun, t=t[m], c=c[m], success=r.success)


def calibrate_turbulence_classes(t, c, geom, ws, p, free=("K0", "tau"), K0=1e-4, tau=10.0, Kmin=5e-8,
                                 t_start=0.0, t_end=None, N=250, bounds=None, shape=None):
    """As calibrate_turbulence but the standard is a size RANGE: class velocities ws and fractions p are
    fixed (known size range and density); K0/tau/Kmin in `free` are fitted.
    With `shape` (from adv_shape) the turbulence is K(t) = Kmin + K0*shape(t); then `free` may contain only
    'K0' and 'Kmin' (tau does not exist)."""
    if shape is not None and "tau" in free:
        raise ValueError("with an ADV shape the free parameters can only be 'K0' and/or 'Kmin'")
    t, c = np.asarray(t, float), np.asarray(c, float)
    m = (t >= t_start) & ((t <= t_end) if t_end is not None else np.ones_like(t, bool))
    tt, cc = t[m], c[m]
    start = dict(K0=K0, tau=tau, Kmin=Kmin)
    bnd = dict(K0=(1e-7, 1e-1), tau=(0.5, 900.0), Kmin=(5e-9, 1e-2))
    if bounds:
        bnd.update(bounds)
    free = list(free)
    def unpack(x):
        q = dict(start); q.update({k: float(np.exp(v)) for k, v in zip(free, x)}); return q
    def resid(x):
        q = unpack(x)
        Kf = shape_K(shape, q["K0"], q["Kmin"]) if shape is not None else None
        return model_c_classes(tt, ws, p, q["K0"], q["tau"], q["Kmin"], geom, N, Kfun=Kf) - cc
    r = least_squares(resid, [np.log(start[k]) for k in free],
                      bounds=([np.log(bnd[k][0]) for k in free], [np.log(bnd[k][1]) for k in free]), diff_step=DIFF_STEP)
    rho1 = np.corrcoef(r.fun[:-1], r.fun[1:])[0, 1]
    infl = (1 + max(rho1, 0)) / (1 - min(max(rho1, 0), 0.95))
    cov = np.linalg.pinv(r.jac.T @ r.jac) * np.sum(r.fun ** 2) / max(len(cc) - len(free), 1) * infl
    out = dict(unpack(r.x), rms=float(np.sqrt(np.mean(r.fun ** 2))), rho1=float(rho1), free=free,
               resid=r.fun, t=tt, c=cc)
    out["sigma_ln"] = dict(zip(free, np.sqrt(np.diag(cov))))
    return out


# ----------------------------------------------------------------------------- particle size distribution (PSD) runs
def _d50(d_um, p):
    """Median of the distribution (interpolated on the cumulative mass fraction at bin mid-points)."""
    return float(np.interp(0.5, np.cumsum(p) - 0.5 * p, d_um))


def _sigma(r, nfree):
    """Lag-1 residual autocorrelation and the autocorrelation-inflated fractional 1-sigma of each free parameter
    (in the fitted, log, variables) from a least_squares result."""
    rho1 = float(np.corrcoef(r.fun[:-1], r.fun[1:])[0, 1])
    infl = (1 + max(rho1, 0)) / (1 - min(max(rho1, 0), 0.95))
    cov = np.linalg.pinv(r.jac.T @ r.jac) * np.sum(r.fun ** 2) / max(len(r.fun) - nfree, 1) * infl
    return rho1, np.sqrt(np.diag(cov))


def load_psd(path, d_col=None, p_col=None, p_type="fraction", d_units="um", d_min=None, d_max=None):
    """Read a PSD csv. One row per size bin; the diameter column and the amount column are `d_col` and `p_col`
    (names or 0-based positions; default: the first and second columns).
    p_type: 'fraction'   amount in each bin (any units: %, volume or mass; rescaled to sum 1),
            'cumulative' cumulative percent/fraction FINER than the diameter (differenced to bin amounts).
    d_units: 'um', 'mm' or 'm' (converted to um).  d_min/d_max keep only that diameter range (um).
    Returns (d_um, p) sorted by diameter, p summing to 1, zero bins removed."""
    df = pd.read_csv(path)
    col = lambda c, k: df.iloc[:, k] if c is None else (df.iloc[:, c] if isinstance(c, int) else df[c])
    d = col(d_col, 0).to_numpy(float) * {"um": 1.0, "mm": 1e3, "m": 1e6}[d_units]
    a = col(p_col, 1).to_numpy(float)
    ok = np.isfinite(d) & np.isfinite(a)
    d, a = d[ok], a[ok]
    o = np.argsort(d); d, a = d[o], a[o]
    if p_type == "cumulative":
        a = np.diff(np.r_[0.0, a])
    elif p_type != "fraction":
        raise ValueError("p_type must be 'fraction' or 'cumulative'")
    keep = a > 0
    if d_min is not None:
        keep &= d >= d_min
    if d_max is not None:
        keep &= d <= d_max
    d, a = d[keep], a[keep]
    if d.size == 0:
        raise ValueError("no PSD bins left after reading / filtering")
    return d, a / a.sum()


def rebin_psd(d_um, p, n_max=12):
    """Merge a fine PSD into at most n_max classes with edges equally spaced in log d; class diameter =
    amount-weighted mean, class fraction = sum.  Fewer classes = faster fits (about 0.3 s per 7-class run)."""
    d_um, p = np.asarray(d_um, float), np.asarray(p, float)
    if d_um.size <= n_max:
        return d_um, p
    edges = np.linspace(np.log(d_um.min()), np.log(d_um.max()) + 1e-9, n_max + 1)
    idx = np.clip(np.digitize(np.log(d_um), edges) - 1, 0, n_max - 1)
    dd, pp = [], []
    for k in range(n_max):
        m = idx == k
        if m.any():
            dd.append(np.sum(d_um[m] * p[m]) / p[m].sum()); pp.append(p[m].sum())
    return np.array(dd), np.array(pp)


def signal_weights(d_um, p, gain_exp=0.0):
    """OBS signal fraction of each class: p_i (d_i/d50)^(-gain_exp), renormalised.  gain_exp = 0 treats the PSD
    amounts as signal fractions (equal sensor response per unit amount); ~1 is typical of backscatter per unit
    mass for particles larger than the wavelength."""
    d_um, p = np.asarray(d_um, float), np.asarray(p, float)
    w = p * (d_um / _d50(d_um, p)) ** (-gain_exp)
    return w / w.sum()


def fit_psd(t, c, geom, d_um, p, which, T, S, coef, K0, tau=30.0, Kmin=5e-8, shape=None, nf=3.0, dp_um=8.0,
            rhos=2650.0, gain_exp=0.0, t_start=0.0, t_end=None, N=250):
    """Fit the settling of a size distribution to the OBS decay.
    which = 'rhos': fit the particle density (solid particles, `nf` fixed, default 3);
    which = 'nf'  : fit the fractal dimension (density `rhos`, primary size `dp_um` fixed).
    Every class settles with the velocity the settling law gives for its own size, so one parameter sets the
    whole velocity spectrum.  Turbulence: K(t) = Kmin + K0*shape(t) if `shape` (from adv_shape) is given, else
    Kmin + K0 exp(-t/tau).  Returns the fitted parameter and its 1-sigma, the class velocities, the velocity at
    d50 (with the range implied by the 1-sigma), the mass-weighted mean velocity, rms, lag-1 autocorrelation and
    residuals."""
    t, c = np.asarray(t, float), np.asarray(c, float)
    m = (t >= t_start) & ((t <= t_end) if t_end is not None else np.ones_like(t, bool))
    tt, cc = t[m], c[m]
    d_um, p = np.asarray(d_um, float), np.asarray(p, float)
    pw = signal_weights(d_um, p, gain_exp)
    d50 = _d50(d_um, p)
    rho_w = density(T, S)
    Kf = shape_K(shape, K0, Kmin) if shape is not None else None
    if which == "rhos":
        to_par = lambda x: rho_w + np.exp(x)
        from_par = lambda v: np.log(v - rho_w)
        lo, hi = [np.log(1.0)], [np.log(6000.0)]
        wvec = lambda v: class_ws(d_um, None, T, S, coef, rhos=v, nf=nf, dp_um=dp_um)
        w50f = lambda v: float(settlingvel(d50, dp_um, nf, coef, v, T, S))
    elif which == "nf":
        to_par = lambda x: float(x)
        from_par = lambda v: v
        lo, hi = [1.2], [3.0]
        wvec = lambda v: class_ws(d_um, None, T, S, coef, rhos=rhos, nf=v, dp_um=dp_um)
        w50f = lambda v: float(settlingvel(d50, dp_um, v, coef, rhos, T, S))
    else:
        raise ValueError("which must be 'rhos' or 'nf'")
    t50 = t50_estimate(tt, cc)
    start = np.nan
    if np.isfinite(t50):
        w0 = geom.fall / t50
        start = invert_ws(w0, d50, which, T, S, coef, rhos=rhos, nf=nf, dp_um=dp_um)
    if not np.isfinite(start):
        start = 1800.0 if which == "rhos" else 2.5
    x0 = [float(np.clip(from_par(start), lo[0] + 1e-6, hi[0] - 1e-6))]
    resid = lambda x: model_c_classes(tt, wvec(to_par(x[0])), pw, K0, tau, Kmin, geom, N, Kfun=Kf) - cc
    r = least_squares(resid, x0, bounds=(lo, hi), diff_step=DIFF_STEP)
    val = to_par(r.x[0])
    rho1, sig = _sigma(r, 1)
    if which == "rhos":
        par_lo, par_hi = rho_w + np.exp(r.x[0] - sig[0]), rho_w + np.exp(r.x[0] + sig[0])
        par_sigma = (par_hi - par_lo) / 2
    else:
        par_sigma = abs(r.x[0]) * sig[0]
        par_lo, par_hi = val - par_sigma, val + par_sigma
    wc = wvec(val)
    w50_lo, w50_hi = sorted([w50f(float(np.clip(par_lo, lo[0] if which == "nf" else rho_w + 1.0, 6000.0 if which == "rhos" else 3.0))),
                             w50f(float(np.clip(par_hi, lo[0] if which == "nf" else rho_w + 1.0, 6000.0 if which == "rhos" else 3.0)))])
    return dict(which=which, value=float(val), sigma=float(par_sigma), ws=wc, d_um=d_um, p=p, d50=d50, w50=w50f(val),
                w50_range=(w50_lo, w50_hi), w_mean=float(np.sum(p * wc)), rms=float(np.sqrt(np.mean(r.fun ** 2))),
                rho1=rho1, resid=r.fun, t=tt, c=cc, success=r.success, weights=pw)


def model_c_psd(t, fit, geom, K0, tau, Kmin, shape=None, N=250):
    """Model curve at times t for the classes and velocities of a fit_psd() / fit_ws_classes() result."""
    Kf = shape_K(shape, K0, Kmin) if shape is not None else None
    return model_c_classes(t, fit["ws"], fit.get("weights", fit["p"]), K0, tau, Kmin, geom, N, Kfun=Kf)


def sensitivity_psd(t, c, geom, d_um, p, which, T, S, coef, K0_list, tau=30.0, Kmin=5e-8, shape=None, **kw):
    """Refit for each K0 in K0_list (same shape / tau / Kmin); shows how much the turbulence amplitude moves the answer."""
    rows = []
    for K0 in K0_list:
        f = fit_psd(t, c, geom, d_um, p, which, T, S, coef, K0, tau, Kmin, shape=shape, **kw)
        rows.append(dict(K0=K0, value=f["value"], w50_mm_s=f["w50"] * 1e3, rms=f["rms"]))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------- calibration file (step 2 -> step 3 handoff)
CAL_FORMAT = "settling_calibration"
CAL_VERSION = 1


def calibration_sensitivity(t, c, geom, d_lo, d_hi, n_classes, rhos_list, T, S, coef, shape=None,
                            free=("K0", "Kmin"), K0=1e-4, tau=15.0, Kmin=5e-6, weights=None, N=200):
    """Recalibrate the turbulence for a range of assumed densities of the standard (the standard's density is
    the main uncertainty of K0).  Returns a DataFrame: rho_standard, K0, tau, Kmin, rms."""
    d, p = size_classes(d_lo, d_hi, n_classes, weights)
    rows, start = [], dict(K0=K0, tau=tau, Kmin=Kmin)
    for rho in rhos_list:
        ws = class_ws(d, None, T, S, coef, rhos=rho)
        r = calibrate_turbulence_classes(t, c, geom, ws, p, free=free, shape=shape, N=N, **start)
        rows.append(dict(rho_standard=rho, K0=r["K0"], tau=r["tau"], Kmin=r["Kmin"], rms=r["rms"]))
    return pd.DataFrame(rows)


def _clean(o):
    if isinstance(o, dict):
        return {k: _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, pd.DataFrame):
        return _clean(o.to_dict(orient="records"))
    return o


def save_calibration(path, cal, geom, shape=None, rpm=None, standard=None, sensitivity=None, files=None, notes=""):
    """Write the turbulence calibration to a json file for step 3 (SettlingFit).
    cal: result of calibrate_turbulence_classes; shape: adv_shape dict (None for the exponential model);
    rpm: mixer rotation rate; standard: dict describing the standard (size range, density, T, S, coef ...);
    sensitivity: DataFrame from calibration_sensitivity; files: dict of source file names; notes: free text."""
    import os
    d = dict(format=CAL_FORMAT, version=CAL_VERSION, created=datetime.datetime.now().isoformat(timespec="seconds"),
             rpm=rpm, geometry=dict(h=geom.h, z_obs=geom.z_obs),
             turb_mode="adv" if shape is not None else "exp",
             K0=cal["K0"], Kmin=cal["Kmin"], tau=None if shape is not None else cal["tau"],
             free=list(cal["free"]), rms=cal["rms"], rho1=cal["rho1"], sigma_ln=cal.get("sigma_ln"),
             shape=(None if shape is None else dict(t=shape["t"], shape=shape["shape"], noise=shape["noise"])),
             standard=standard, K0_vs_standard_density=sensitivity, files=files, notes=notes)
    folder = os.path.dirname(os.path.abspath(path))
    os.makedirs(folder, exist_ok=True)
    with open(path, "w") as f:
        json.dump(_clean(d), f, indent=1)
    return path


def load_calibration(path):
    """Read a calibration file. Returns a dict with K0, Kmin, tau (30 s placeholder in 'adv' mode, where it is
    unused), shape (adv_shape-style dict or None), geom (Geometry), rpm, turb_mode and the stored diagnostics."""
    with open(path) as f:
        d = json.load(f)
    if d.get("format") != CAL_FORMAT:
        raise ValueError(f"{path} is not a settling calibration file")
    sh = d.get("shape")
    d["shape"] = None if sh is None else dict(t=np.asarray(sh["t"], float), shape=np.asarray(sh["shape"], float), noise=sh["noise"])
    d["geom"] = Geometry(h=d["geometry"]["h"], z_obs=d["geometry"]["z_obs"])
    d["tau_used"] = 30.0 if d.get("tau") is None else d["tau"]
    return d


def check_calibration(cal, geom=None, rpm=None, tol=0.02):
    """List of warnings if the calibration does not apply to the stated setup (empty list = consistent)."""
    w = []
    if geom is not None:
        g = cal["geom"]
        if abs(geom.h / g.h - 1) > tol:
            w.append(f"water depth {geom.h:.4f} m differs from the calibration ({g.h:.4f} m)")
        if abs(geom.z_obs / g.z_obs - 1) > tol:
            w.append(f"sensor height {geom.z_obs:.4f} m differs from the calibration ({g.z_obs:.4f} m)")
    if rpm is not None and cal.get("rpm") is not None and rpm != cal["rpm"]:
        w.append(f"mixer rate {rpm} rpm differs from the calibration ({cal['rpm']} rpm)")
    return w


def predict_standard(t, geom, cal, d_lo, d_hi, rho, T, S, coef, n_classes=7, weights=None, N=250):
    """Model decay of a standard (size range, density) with a stored calibration and NO fitting; compare with data
    to test whether the calibration transfers to another standard."""
    d, p = size_classes(d_lo, d_hi, n_classes, weights)
    ws = class_ws(d, None, T, S, coef, rhos=rho)
    return model_c_classes(t, ws, p, cal["K0"], cal["tau_used"], cal["Kmin"], geom, N,
                           Kfun=shape_K(cal["shape"], cal["K0"], cal["Kmin"]) if cal["shape"] is not None else None)


# ----------------------------------------------------------------------------- one run, start to finish (step 3 helper)
def fit_run(obs_file, psd_file, geom, cal, T, S, coef, what="rhos", nf=3.0, dp_um=8.0, rhos_for_nf=2650.0, beta=1.3,
            gain_exp=0.0, psd_kw=None, n_max_classes=12, t_start=0.0, t_end=None, K0_factor=1.4, do_sensitivity=True,
            N=250, name=None, load_kw=None):
    """Fit one OBS decay with one PSD using a stored calibration (from load_calibration).
    what: 'rhos' (density, solid particles), 'nf' (fractal dimension) or 'ws50' (velocity at the median size, w ~ d^beta).
    If do_sensitivity, the fit is repeated with K0 divided and multiplied by K0_factor.
    Returns (row, fit, data): row = dict for a results table, fit = the fit dict, data = dict(t, c, d_um, p, name)."""
    d = load_decay(obs_file, **(load_kw or {}))
    t, c = d["t"], d["c"]
    dd, pp = load_psd(psd_file, **(psd_kw or {}))
    dd, pp = rebin_psd(dd, pp, n_max_classes)
    K0, tau, Kmin, shape = cal["K0"], cal["tau_used"], cal["Kmin"], cal["shape"]

    def one(K0_):
        if what in ("rhos", "nf"):
            return fit_psd(t, c, geom, dd, pp, what, T, S, coef, K0_, tau, Kmin, shape=shape, nf=nf, dp_um=dp_um,
                           rhos=rhos_for_nf, gain_exp=gain_exp, t_start=t_start, t_end=t_end, N=N)
        f = fit_ws_classes(t, c, geom, dd, pp, K0_, tau, Kmin, beta=beta, shape=shape, t_start=t_start, t_end=t_end, N=N)
        f.update(value=f["ws50"], sigma=f["ws50_sigma"], w50=f["ws50"], weights=pp, d_um=dd, w_mean=float(np.sum(pp * f["ws"])),
                 w50_range=(f["ws50"] - f["ws50_sigma"], f["ws50"] + f["ws50_sigma"]))
        return f

    fit = one(K0)
    scale = 1e3 if what == "ws50" else 1.0
    row = dict(run=name or str(obs_file), what=what,
               units={"rhos": "kg/m3", "nf": "-", "ws50": "mm/s"}[what],
               value=fit["value"] * scale, sigma=fit["sigma"] * scale,
               w50_mm_s=fit["w50"] * 1e3, w50_lo_mm_s=fit["w50_range"][0] * 1e3, w50_hi_mm_s=fit["w50_range"][1] * 1e3,
               w_mean_mm_s=fit["w_mean"] * 1e3, d50_um=fit["d50"], rms=fit["rms"], rho1=fit["rho1"], K0=K0, Kmin=Kmin)
    if do_sensitivity:
        lo, hi = one(K0 / K0_factor), one(K0 * K0_factor)
        row.update(value_K0_low=lo["value"] * scale, value_K0_high=hi["value"] * scale, K0_factor=K0_factor)
    return row, fit, dict(t=t, c=c, d_um=dd, p=pp, name=row["run"])
