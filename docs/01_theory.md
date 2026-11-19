# 1. Theory of the fit

This page explains what the model assumes, how a settling velocity is extracted from a concentration decay, and what limits the answer. For how to run it, see [02_user_guide.md](02_user_guide.md); for every option, see [03_reference.md](03_reference.md).

## 1.1 The measurement and the idea

A mixer keeps a suspension fully mixed in a tank. An optical backscatter sensor (OBS) at height $z_{obs}$ above the bed (depth $z_s = h - z_{obs}$ below the free surface, $h$ = water depth) records the total concentration. When the mixer stops (t = 0) particles settle out of the volume above the sensor and the signal falls from its initial value $C_o$ to a background value $C_{bg}$.

If the water were perfectly still, a particle of settling velocity $w$ would leave the sensor level at $t = (h - z_{obs})/w$, and the decay would be a step (one size) or a smoothed ramp (a size distribution). In practice the tank is not still: turbulence from the mixer decays over roughly 10–35 s and keeps mixing particles vertically while they settle. The decay curve therefore contains both settling and residual mixing. The model separates the two by describing the mixing explicitly with a diffusivity $K(t)$ that is calibrated on a **size standard** (known size and density), so that the settling velocity of an unknown population can then be found from its decay alone.

Typical numbers (75–90 µm standard, $\rho_s$ = 1200 kg m⁻³, $h$ = 0.178 m, sensor 0.064 m above the bed): $w \approx 0.67$ mm s⁻¹, so clearing the sensor takes about 170 s; turbulence is gone after about 30 s; $K_0 \approx 10^{-4}$ m² s⁻¹. The diffusive spreading during the turbulent period, $\sqrt{4\int K\,dt}$ with $\int K\,dt \approx K_0\tau \approx 1.4\times10^{-3}$ m², is about 7 cm, comparable to the 11 cm between the sensor and the free surface. Free-space (erfc) solutions that ignore the surface and bed are therefore not accurate here, which is why the model includes the tank boundaries.

## 1.2 Governing equation

With $z$ measured **downward from the free surface** ($z = 0$ surface, $z = h$ bed), the concentration $c(z,t)$ of a size class with settling velocity $w$ obeys

$$\frac{\partial c}{\partial t} + \frac{\partial}{\partial z}\left(w\,c - K(t)\,\frac{\partial c}{\partial z}\right) = 0 .$$

* $K$ is uniform in depth and depends only on time after the mixer stops.
* **Free surface** ($z = 0$): no particle flux, $w c - K\,\partial_z c = 0$.
* **Bed** ($z = h$): particles that reach the bed stay there, flux $= D\,w\,c$ with $D = 1$ (perfect deposition; $D = 0$ would be a reflecting bed). The bed absorbs settled mass and the model tracks it, so mass is conserved.
* The **initial state** is the steady profile in a stirred tank with diffusivity $K_0 = K(0)$ and zero net flux:

$$c(z,0) \;\propto\; \exp\!\left[\frac{w\,(z-h)}{K_0}\right],$$

normalised to the depth-mean concentration of the class. The Péclet number $Pe = w h/K_0$ says how uniform this is: close to uniform when $Pe \ll 1$, bottom-heavy when $Pe \gtrsim 1$. For the standard above, $Pe \approx 1.2$, so the stirred profile is measurably non-uniform and the initial value at the sensor differs from the depth mean. The model handles this because it normalises by its own initial value at the sensor (section 1.6).

## 1.3 The turbulence model $K(t)$

Two forms are available.

**Exponential:** $K(t) = K_{min} + K_0\,e^{-t/\tau}$, with $K_0$, $\tau$ and $K_{min}$ fitted on the standard.

**ADV shape:** $K(t) = K_{min} + K_0\,s(t)$, where $s(t)$ is a dimensionless decay shape measured with an acoustic Doppler velocimeter (ADV) in the same tank at the same mixer rate, with $s(0) = 1$. This removes the strong trade-off between $K_0$ and $\tau$ that a decay curve alone cannot resolve.

**How the shape is obtained** (`adv_shape`):

1. The ADV decay record is split into 2 s windows (50 % overlap).
2. In each window the linear trend of each velocity component is removed, so the slowly varying mean flow and instrument offsets are excluded, and the fluctuation variance $E = \overline{v_1'^2 + v_2'^2 + v_3'^2}$ is computed.
3. The noise floor $E_{noise}$ (the median of $E$ after 60 s, when the water is effectively still and only instrument noise remains) is subtracted.
4. The diffusivity is assumed proportional to the rms fluctuation velocity, $K \propto u' \propto \sqrt{E - E_{noise}}$ (an eddy diffusivity is a velocity scale times a length scale; the length scale is taken as constant in time). A short running mean is applied.
5. The shape is set to zero after `t_cut_s` (30 s by default), where the ADV signal is indistinguishable from noise, and normalised to 1 in the first window.

Only the **shape** comes from the ADV. The amplitude $K_0$ cannot be obtained from point velocities without a mixing length, so it is calibrated (section 1.8). $K_{min}$ represents a persistent weak mixing that the ADV cannot resolve (a floor of a few $10^{-6}$ m² s⁻¹ corresponds to velocity fluctuations of about 0.1 mm s⁻¹ at a scale of a few centimetres, far below the ADV's noise level of roughly 7–12 mm s⁻¹); it is an effective model parameter, not a measured quantity. If $K \propto u'$ is wrong, for example if the mixing length changes as the turbulence decays, the shape is biased; the exponential option is the alternative that avoids this assumption but loses the ADV constraint.

## 1.4 Numerical solution

The domain $0 \le z \le h$ is divided into $N$ equal cells (200–600 depending on the function). Advective fluxes at cell faces use a second-order Van Leer flux limiter (monotone, with little numerical diffusion, which matters late in the decay when $K$ is small and the settling front is nearly sharp); diffusive fluxes are central differences. The semi-discrete system is integrated in time with the implicit variable-step BDF method (`scipy.integrate.solve_ivp`). For several size classes the classes are integrated together in one system (they do not interact), with a coloured finite-difference Jacobian that needs five perturbed evaluations for the whole matrix. The multi-class solution reproduces the sum of single-class solutions to about $10^{-7}$. A first-order exponential-fitting (Scharfetter–Gummel) scheme is also available in `solve_single(..., scheme="sg")` for comparison. The solution should be insensitive to $N$; compare $N$ = 200 and 400 if in doubt.

## 1.5 Settling velocity and the settling law

A single parameter sets the velocity of every size class through a settling law (Ferguson & Church, extended to fractal aggregates):

$$w = \frac{g\,R\,d^{\,n_f-1}}{C_1\,\nu\,d_p^{\,n_f-3} + C_2\sqrt{g\,R\,d^{\,n_f}\,d_p^{\,n_f-3}}},\qquad R = \frac{\rho_s-\rho_w}{\rho_w}.$$

* $d$ is the particle (floc) diameter, $d_p$ the primary particle diameter, $n_f$ the fractal dimension ($n_f = 3$ is a solid particle, for which $d_p$ drops out).
* $\rho_w(T, S)$ and the kinematic viscosity $\nu(T)$ come from the water temperature and salinity.
* $(C_1, C_2)$ are set by the coefficient choice: 1 flocs matching Stokes at the primary particles (20, 1.258), 2 viscous flocs (120, 0), 3 sieve-sized sand or smooth spheres (18, 0.866), 4 nominal-diameter natural sand (20, 0.91), 5 sphere (18, 0.548).

A fitted $w$ is converted back to a density (solid particles, $n_f$ fixed) or to a fractal dimension (flocs, density and primary size fixed) by solving this relation. The inversion returns "no solution" if $w$ lies outside what the law can produce for the given size.

## 1.6 What the sensor sees: observation model and size distributions

The model returns the concentration in the cell containing the sensor depth $z_s$. A Gaussian sensor footprint is available in the solver (`sigma_s`) but is off by default.

Data and model are both normalised so their absolute scale does not matter:

* data: $c(t) = (\mathrm{NTU}(t) - C_{bg})/(C_o - C_{bg})$, using the initial plateau $C_o$ and the background $C_{bg}$ identified in step 1;
* model: divided by its own value at $t = 0$.

For a size distribution the classes are independent (the equation is linear), so the signal is the sum over classes $i$ of the class concentration at the sensor, each starting from its own stirred profile with depth-mean $p_i$:

$$c_{obs}(t) = \sum_i c_i(z_s, t).$$

$p_i$ is the fraction of the initial OBS signal carried by class $i$. If the PSD file gives mass or volume fractions and the sensor responds differently per unit mass at different sizes (backscatter per unit mass falls with size for particles larger than the wavelength), a response weighting $p_i \propto \text{amount}_i\,(d_i/d_{50})^{-\mathrm{gain\_exp}}$ can be applied. The default, `gain_exp = 0`, treats the PSD amounts as signal fractions. A single size is a PSD with one class.

## 1.7 Estimation

The unknown (density, fractal dimension, or a velocity scale) is found by non-linear least squares on the normalised decay, with parameters fitted in logarithmic form where they must stay positive. Finite-difference derivatives use a relative step of $10^{-3}$ in the log-parameters, because the ODE solver output carries about $10^{-5}$ integration noise and the default step would give unreliable uncertainties.

**Uncertainty.** The quoted 1-σ comes from the Jacobian and the residual variance, inflated by $(1+\rho_1)/(1-\rho_1)$ where $\rho_1$ is the lag-1 autocorrelation of the residuals ($\rho_1$ clipped to 0–0.95). Residuals in these decays are autocorrelated ($\rho_1 \approx 0.4$ for the standard), so uncorrelated-noise errors would be too small. The 1-σ is a **fit uncertainty only**. It does not include errors in $K_0$ and $K_{min}$, the geometry, the PSD or the settling law.

**Three things can be fitted** to a run:

| Mode | Unknown | How class velocities follow |
|---|---|---|
| `rhos` | density of solid particles ($n_f$ fixed) | settling law at each class diameter |
| `nf` | fractal dimension of flocs (density and $d_p$ fixed) | settling law at each class diameter |
| `ws50` | velocity at the median size | $w_i = w_{50}(d_i/d_{50})^{\beta}$, no density or $n_f$ assumption |

`ws50` is the least assumption-laden summary of the settling (it needs only the slope $\beta$ of the velocity spectrum: about 2 in the Stokes range, about 1.3 for 75–90 µm sand-like spheres in this law). `rhos` and `nf` add the settling law.

## 1.8 Why and how the turbulence is calibrated

On a decay curve alone, a faster-settling population with more mixing and a slower population with less mixing can look alike, so $w$ and $K$ are partly interchangeable. The calibration breaks this by using a standard whose settling velocities are known from its size range and density. With those velocities **fixed**, only the turbulence parameters ($K_0$ and $K_{min}$, or $K_0$, $\tau$ and $K_{min}$) are fitted. The resulting $K(t)$ is then held fixed when fitting unknown populations in the same tank at the same mixer rate.

The calibration is exactly as good as the standard's description. $K_0$ trades off strongly against the assumed density: for the 75–90 µm standard, changing the density by ±50 kg m⁻³ changes $K_0$ by a factor of about 2.5 (4.2×10⁻⁵ to 2.6×10⁻⁴ m² s⁻¹) while the fit quality barely changes. The calibration notebook tabulates this, and the fit notebook propagates a $K_0$ uncertainty factor (`K0_factor`) into every result so that it is visible next to the fit's own ±.

The real test of a calibration is to use it to **predict** the decay of a second standard of different size or density with no refitting. A prediction rms near the noise level shows that the turbulence model transfers between particle sizes; a clearly larger rms points to the assumed density or size distribution of one of the standards, or to the 1-D model.

## 1.9 Assumptions and limitations

* One-dimensional, with depth-uniform diffusivity. The tank's spin-down circulation (secondary flows after the mixer stops) is not represented.
* $K(t)$ follows an exponential or an ADV-derived shape with $K \propto u'$; the ADV is a point measurement, not a tank average.
* Particles do not interact, flocculate or break up during the decay, and the OBS responds linearly. Size-dependent sensor response is only handled through `gain_exp`.
* Bed: perfect deposition; no resuspension after the mixer stops.
* The calibration applies to one mixer rate in one tank setup (water depth, sensor position).
* For the 75–90 µm standard the best fit still leaves a residual rms of about 0.08 against sensor noise of about 0.04, with lag-1 autocorrelation of about 0.4. Something systematic remains (the density or size distribution assumed for the standard, or the 1-D model); do not read the quoted ± as the full uncertainty.
