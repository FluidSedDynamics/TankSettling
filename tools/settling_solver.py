"""
1D settling in a mixing tank after the mixer stops.

    dc/dt + d/dz( w c - K(t) dc/dz ) = 0 ,   z downward from the free surface (0) to the bed (h)

* Finite volume, exponential-fitting (Scharfetter-Gummel) face flux: exact for steady
  advection-diffusion between cell centres, so it stays stable and monotone for any
  cell Peclet number (important once K(t) has decayed and the front is nearly sharp).
* Time integration: scipy BDF on the semi-discrete linear system with tridiagonal Jacobian.
* Initial state: the stirred steady profile (zero flux at surface and bed, uniform K0):
  c ~ exp(w (z - h) / K0), i.e. uniform when Pe = w h / K0 << 1, bottom-heavy otherwise.
* Surface BC: zero flux.  Bed BC after mixer-off: F = D * w * c_bed  (D=1 absorbing / perfect
  deposition, D=0 reflecting).
* Observation: Gaussian-weighted average about depth z_s (sensor footprint sigma_s), or a
  single cell if sigma_s = 0.
* Multi-fraction: the problem is linear and fractions do not interact, so the total signal is
  sum_i k_i * p_i * c_obs(t; w_i).
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import diags, csr_matrix


# ----------------------------------------------------------------------------- turbulence
def K_exp(K0, tau, Kmin=0.0):
    """Uniform-in-z diffusivity decaying exponentially after mixer-off."""
    return lambda t: K0 * np.exp(-t / tau) + Kmin


def theta_exp(K0, tau, Kmin=0.0):
    """theta(t) = integral_0^t K dt'  (the 'diffusive time'); front width = sqrt(4 theta)."""
    return lambda t: K0 * tau * (1.0 - np.exp(-t / tau)) + Kmin * t


# ----------------------------------------------------------------------------- discretisation
def _grid(h, N):
    dz = h / N
    zc = (np.arange(N) + 0.5) * dz
    return dz, zc


def _bernoulli_flux_coeffs(w, K, dz):
    """Face flux F = a*c_i - b*c_{i+1} for w>0 (downward), exact for constant F across the face.
    F = w (c_i - c_{i+1} e^{-a}) / (1 - e^{-a}),  a = w dz / K."""
    K = max(K, 1e-300)
    a = w * dz / K
    if a > 700.0:                      # pure upwind
        return w, 0.0
    if a < 1e-8:                       # pure diffusion limit
        return K / dz + 0.5 * w, K / dz - 0.5 * w
    e = np.exp(-a)
    den = -np.expm1(-a)                # 1 - e^{-a}
    return w / den, w * e / den


def _rhs_matrix(w, K, dz, N, D):
    """dc/dt = A c with A tridiagonal. Face fluxes (downward positive):
       F_0 = 0, F_{i+1/2} = alpha c_i - beta c_{i+1}, F_N = D w c_{N-1}."""
    alpha, beta = _bernoulli_flux_coeffs(w, K, dz)
    # interior faces i+1/2 for i=0..N-2
    # dc_i/dt = -(F_{i+1/2} - F_{i-1/2}) / dz
    main = np.zeros(N)
    lower = np.zeros(N - 1)   # coefficient of c_{i-1} in row i  (i=1..N-1)
    upper = np.zeros(N - 1)   # coefficient of c_{i+1} in row i  (i=0..N-2)
    # outflow from cell i through its bottom face: alpha*c_i - beta*c_{i+1}
    main[:-1] += -alpha / dz
    upper[:] += beta / dz
    # inflow to cell i through its top face (i>=1): alpha*c_{i-1} - beta*c_i
    lower[:] += alpha / dz
    main[1:] += -beta / dz
    # bed face
    main[-1] += -D * w / dz
    return diags([lower, main, upper], [-1, 0, 1], format="csr")


# ----------------------------------------------------------------------------- initial state
def stirred_profile(w, K0, h, N, mean=1.0):
    """Steady zero-flux profile under uniform K0, normalised to depth-mean `mean`."""
    dz, zc = _grid(h, N)
    if w <= 0:
        return np.full(N, mean)
    c = np.exp(w * (zc - h) / K0)
    return c * mean / c.mean()


# ----------------------------------------------------------------------------- observation
def sensor_weights(zc, z_s, sigma_s, dz):
    if sigma_s <= 0:
        wts = np.zeros_like(zc)
        wts[np.argmin(np.abs(zc - z_s))] = 1.0
        return wts
    wts = np.exp(-0.5 * ((zc - z_s) / sigma_s) ** 2)
    return wts / wts.sum()



# ----------------------------------------------------------------------------- limited (2nd order) scheme
def _rhs_limited(c, w, K, dz, D):
    """Van Leer-limited advective face values (w>0, downward) + central diffusion.
    Returns dc/dt and the bed deposition flux.  Much less numerical diffusion than the
    exponential-fitting (upwind-like) flux when K is small (late times)."""
    N = c.size
    dc = c[1:] - c[:-1]
    num = np.empty(N - 1)
    num[0] = 0.0
    num[1:] = c[1:-1] - c[:-2]
    lim = (num * np.abs(dc) + np.abs(num) * dc) / (np.abs(dc) + np.abs(num) + 1e-300)
    F = np.empty(N + 1)
    F[0] = 0.0
    F[1:N] = w * (c[:-1] + 0.5 * lim) - K * dc / dz
    F[N] = D * w * c[-1]
    return -(F[1:] - F[:-1]) / dz, F[N]


def _solve_limited(w, Kfun, h, z_s, t_eval, N, D, K0_init, sigma_s, c_init, rtol, atol, return_profiles):
    dz, zc = _grid(h, N)
    c0 = stirred_profile(w, K0_init, h, N) if c_init is None else np.asarray(c_init, float)
    wts = sensor_weights(zc, z_s, sigma_s, dz)

    def f(t, y):
        dcdt, dep = _rhs_limited(y[:-1], w, Kfun(t), dz, D)
        return np.concatenate([dcdt, [dep]])

    pat = diags([np.ones(N + 1 - abs(k)) for k in range(-2, 3)], list(range(-2, 3)), format="csr")
    sol = solve_ivp(f, (0.0, float(np.max(t_eval))), np.concatenate([c0, [0.0]]), method="BDF",
                    t_eval=t_eval, rtol=rtol, atol=atol, jac_sparsity=pat)
    if not sol.success:
        raise RuntimeError(sol.message)
    C = sol.y[:-1, :]
    c_obs = wts @ C
    if return_profiles:
        return c_obs, (zc, C, sol.y[-1, :] / h)
    return c_obs

# ----------------------------------------------------------------------------- solver
def solve_single(w, Kfun, h, z_s, t_eval, N=400, D=1.0, K0_init=None, sigma_s=0.0,
                 c_init=None, rtol=1e-8, atol=1e-10, return_profiles=False, scheme="limited"):
    """Single size fraction, depth-mean initial concentration 1.

    w        settling velocity (m/s, >0 downward)
    Kfun     callable K(t) (m^2/s), t=0 at mixer-off
    scheme   'limited' (default, 2nd order Van Leer) or 'sg' (exponential fitting, first order)
    K0_init  diffusivity of the stirred state (default Kfun(0))
    Returns c_obs(t_eval) [and (zc, profiles, deposited_mass_fraction(t)) if requested].
    """
    if K0_init is None:
        K0_init = Kfun(0.0)
    if scheme == "limited":
        return _solve_limited(w, Kfun, h, z_s, t_eval, N, D, K0_init, sigma_s, c_init, rtol, atol,
                              return_profiles)
    dz, zc = _grid(h, N)
    c0 = stirred_profile(w, K0_init, h, N) if c_init is None else np.asarray(c_init, float)
    wts = sensor_weights(zc, z_s, sigma_s, dz)

    def f(t, y):
        return _rhs_matrix(w, Kfun(t), dz, N, D) @ y

    def jac(t, y):
        return _rhs_matrix(w, Kfun(t), dz, N, D)

    # augment with deposited mass so the mass balance can be checked
    def f_aug(t, y):
        c = y[:-1]
        return np.concatenate([f(t, c), [D * w * c[-1]]])

    def jac_aug(t, y):
        A = _rhs_matrix(w, Kfun(t), dz, N, D).tolil()
        A.resize((N + 1, N + 1))
        A[N, N - 1] = D * w
        return csr_matrix(A)

    y0 = np.concatenate([c0, [0.0]])
    sol = solve_ivp(f_aug, (0.0, float(np.max(t_eval))), y0, method="BDF", jac=jac_aug,
                    t_eval=t_eval, rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(sol.message)
    C = sol.y[:-1, :]
    c_obs = wts @ C
    if return_profiles:
        dep = sol.y[-1, :] / h   # deposited mass per unit depth-mean (initial mean = 1)
        return c_obs, (zc, C, dep)
    return c_obs


def solve_mixture(ws, p, k, Kfun, h, z_s, t_eval, **kw):
    """Total sensor signal for fractions with settling velocities ws, mass fractions p, OBS gains k."""
    out = np.zeros(len(t_eval))
    for wi, pi, ki in zip(ws, p, k):
        out += ki * pi * solve_single(wi, Kfun, h, z_s, t_eval, **kw)
    return out


# ----------------------------------------------------------------------------- analytic front
def erfc_front(w, theta, z_s, t):
    """Leading-order far-boundary solution c/c0 = 0.5 erfc((w t - z)/sqrt(4 theta)); error O(l_d/z)."""
    from scipy.special import erfc
    return 0.5 * erfc((w * t - z_s) / np.sqrt(4.0 * np.maximum(theta, 1e-300)))


# ----------------------------------------------------------------------------- multi-class (vectorised)
def _rhs_limited_multi(C, w, K, dz, D):
    """C: (M,N) class concentrations, w: (M,1).  Same scheme as _rhs_limited, all classes at once."""
    M, N = C.shape
    dc = C[:, 1:] - C[:, :-1]
    num = np.zeros((M, N - 1))
    num[:, 1:] = C[:, 1:-1] - C[:, :-2]
    lim = (num * np.abs(dc) + np.abs(num) * dc) / (np.abs(dc) + np.abs(num) + 1e-300)
    F = np.empty((M, N + 1))
    F[:, 0] = 0.0
    F[:, 1:N] = w * (C[:, :-1] + 0.5 * lim) - K * dc / dz
    F[:, N] = D * w[:, 0] * C[:, -1]
    return -(F[:, 1:] - F[:, :-1]) / dz, F[:, N]


def solve_classes(ws, p, Kfun, h, z_s, t_eval, N=300, D=1.0, K0_init=None, sigma_s=0.0,
                  rtol=1e-7, atol=1e-9, k=None, return_classes=False):
    """Several size classes (settling velocities ws, signal/mass fractions p) in ONE BDF integration.

    Each class starts from its own stirred steady profile (depth-mean p_i) and obeys the same K(t).
    Returns the sensor signal  sum_i k_i * c_i(z_s, t)  (k_i = 1 by default); the t=0 value is
    sum_i p_i * (stirred-profile value at the sensor)."""
    ws = np.atleast_1d(np.asarray(ws, float)); p = np.asarray(p, float)
    M = ws.size
    if k is None:
        k = np.ones(M)
    if K0_init is None:
        K0_init = Kfun(0.0)
    dz, zc = _grid(h, N)
    wts = sensor_weights(zc, z_s, sigma_s, dz)
    C0 = np.array([stirred_profile(wi, K0_init, h, N, mean=pi) for wi, pi in zip(ws, p)])
    wcol = ws[:, None]
    NM = M * N

    def f(t, y):
        C = y[:NM].reshape(M, N)
        dC, dep = _rhs_limited_multi(C, wcol, Kfun(t), dz, D)
        return np.concatenate([dC.ravel(), dep])

    # coloured finite-difference Jacobian: classes are decoupled and the stencil is +-2 cells, so
    # 5 perturbed evaluations give the whole matrix (scipy's generic sparse FD is far slower).
    NT = NM + M
    colors = []
    for g in range(5):
        cols = np.array([i * N + n for i in range(M) for n in range(g, N, 5)])
        rows, cc = [], []
        for col in cols:
            i, n = divmod(col, N)
            for r in range(max(0, n - 2), min(N, n + 3)):
                rows.append(i * N + r); cc.append(col)
            if n == N - 1:
                rows.append(NM + i); cc.append(col)
        colors.append((cols, np.array(rows), np.array(cc)))

    def jac(t, y):
        f0 = f(t, y)
        data_r, data_c, data_v = [], [], []
        for cols, rows, cc in colors:
            hstep = 1e-7 * np.maximum(np.abs(y[cols]), 1e-3 / M) 
            yp = y.copy(); yp[cols] += hstep
            df = f(t, yp) - f0
            hv = np.zeros(NT); hv[cols] = hstep
            data_r.append(rows); data_c.append(cc); data_v.append(df[rows] / hv[cc])
        return csr_matrix((np.concatenate(data_v), (np.concatenate(data_r), np.concatenate(data_c))),
                          shape=(NT, NT))

    y0 = np.concatenate([C0.ravel(), np.zeros(M)])
    sol = solve_ivp(f, (0.0, float(np.max(t_eval))), y0, method="BDF", t_eval=t_eval,
                    rtol=rtol, atol=atol, jac=jac)
    if not sol.success:
        raise RuntimeError(sol.message)
    Cc = sol.y[:NM, :].reshape(M, N, -1)
    per = np.einsum("n,mnt->mt", wts, Cc)
    sig = (k[:, None] * per).sum(axis=0)
    return (sig, per) if return_classes else sig
