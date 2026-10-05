# The science behind Ink Forest

This document explains what every element of the Ink Forest picture means, for readers with a background in cosmology or statistics. The [README](README.md) gives the short version.

The game is one loop: **observe** a sightline, **infer** the structure between the observed sightlines, **revise** the map and the next choice. Everything on screen is one of three kinds: **measured**, **inferred** or **presentation** (§7).

## 1. The plane and the sightlines

**The plane.**
- The universe in the game is one plane of synthetic Lyα transmitted flux from a PRIYA simulation at redshift z = 3 ([DATA.md](DATA.md)).
- It holds 480 parallel sightlines spaced 0.25 Mpc/h apart.
- Each sightline has 1389 pixels of 10 km/s, spanning the 120 Mpc/h periodic box.

**The survey.**
- 25 candidate sightlines, 1.5 Mpc/h apart, spanning 37.5 Mpc/h across.
- The player may observe 10 of them.

**What is drawn.**
- The sheet shows one 40 Mpc/h stretch of depth, 40–80 Mpc/h along the line of sight, drawn downward.
- Transverse position runs across the sheet.

**Notation.**
- F is the transmitted flux; A = 1 − F is the absorption.
- δ = A − ⟨A⟩ is the absorption contrast. ⟨A⟩ = 0.324 is the mean over sightlines outside the survey.

## 2. Measured: the forest ribbons

**Observation model.**
- An observed sightline is the simulated flux plus Gaussian white noise, σ_F = 0.02 per 20 km/s.
- The noise realisation is fixed per sightline, so a sightline always says the same thing.
- There is no instrumental resolution, continuum error or metal-line contamination.

**The ribbon.**
- Each observed sightline is drawn as a ribbon built only from its own observed spectrum, at full resolution.
- A hairline where the gas is transparent; wider and darker with absorption.
- A broad-absorption term (the line's own absorption smoothed over 1 Mpc/h) widens long saturated stretches, so a damped absorber looks different from a cluster of narrow lines.

The ribbons are the evidence: they are never smoothed toward the inferred field, and nothing else is drawn as sharply.

## 3. Inferred: the sheet between the lines

**Target.** The sheet estimates T, the absorption contrast δ smoothed with a Gaussian of σ = 1.5 Mpc/h in both directions. That scale equals the sightline spacing: the map claims no structure finer than the spacing.

**Estimator.**
- The exact Gaussian conditional mean, E[T | observed sightlines], treating δ as a stationary Gaussian field with white observational noise.
- The observed sightlines are full periodic lines. The covariance is therefore block-circulant along the line of sight, and the solution decouples into one small linear system per line-of-sight Fourier mode (|q| ≤ 60, which carries all of T's power).
- The browser solves this exactly each time the set of observed lines changes, in a few milliseconds.

**Covariance.**
- The two-point function of δ measured from 250 other sightlines of the same plane, all well outside the survey region, so no surveyed sightline informs its own prior.
- Symmetrised.
- Tapered smoothly to zero between 15 and 30 Mpc/h of transverse separation, which suppresses estimation noise at large lags.

**Signed field.**
- What is drawn is s = E[T] / σ_T, in units of the prior standard deviation.
- s > 0: more absorption than average, i.e. denser gas. Drawn as graphite ink.
- s < 0: less absorption than average, i.e. emptier regions. Drawn as lifted, paler paper.
- Both signs are shown with equal care. Their tones are equalised because the field is skewed: walls are deep, empty regions shallow.

**Shrinkage, not hiding.** Far from any observed sightline the conditional mean relaxes toward zero, the cosmic mean, by itself. The renderer never multiplies the field's amplitude by confidence.

## 4. Uncertainty and support

**Support** at each point is c = 1 − Var[T | data] / Var[T]: the fraction of the prior variance explained by the observations.

**Drawn as material certainty, not as a separate overlay.**
- Low support: wet, grainy, soft-edged ink.
- High support: dry, smooth ink with a crisp edge at |s| = 0.75.
- Bare paper: unconstrained.
- There are no uncertainty boundaries, fog or masks.

**In the debrief,** "unresolved" means support below 0.30.

## 5. Damped absorbers and setting a line aside

**What they are.** A damped Lyα absorber (DLA) is a single dense, largely neutral cloud. Its saturated core and damping wings absorb all the light over thousands of km/s of one spectrum, far more than the forest itself.

**Why it misleads.** Read as ordinary forest, that stretch is an enormous absorption excess on one line. The Gaussian estimator propagates it to the neighbouring lines as a broad false wall, and it can break up genuine structure nearby.

**The rule.** Ink Forest uses a fixed rule, applied to the line's own observed spectrum only (pairs of pixels averaged to 20 km/s):

1. A = 1 − (mean flux over ±2 samples).
2. A damped core is a run of A > 0.93 lasting at least 40 samples (800 km/s).
3. The masked stretch is the core extended along the line's own absorption while A > 0.3.

**Setting aside.** Pressing and holding a landed line sets aside its masked stretch.
- The map is recomputed by **exact conditioning on all other pixels**: the masked pixels are given infinite noise, and the result is computed exactly with a low-rank update over all line-of-sight modes.
- The rest of that sightline still informs the map.
- Support drops where information was removed, instead of the gap being filled.
- Holding again restores the ordinary reading.
- On a line where the rule finds no damped stretch, holding does nothing.

**No verdict.** The game never tells the player whether a line contains a DLA.

**Caveat.** The rule sees one spectrum only. A long saturated stretch can also be a genuine dense large-scale structure, which neighbouring sightlines would share. What singles out a lone dense cloud is that its neighbours do not show it, so observing a neighbour is the player's test.

## 6. The ride: what the DIVE shows, and when

The ride is causal:
- Everything drawn or heard at a given moment depends only on the part of the spectrum the drop has already passed.
- Smoothing looks backward only.
- The pace (slower through absorption) is a fixed per-bin rule, not normalised by the whole line.
- An absorber is confirmed, and sounds, only once the drop has left it.

**No early information.**
- The map is not updated until the ride is complete.
- The candidate lights never reveal anything about a sightline before it is chosen: unobserved spectra cannot affect what is on screen.

## 7. Measured, inferred, presentation

| Element | Category |
|---|---|
| Ribbon shape and darkness | **Measured**: the observed spectrum (simulated, with noise) |
| Graphite and lifted washes between the ribbons | **Inferred**: Gaussian conditional mean of the smoothed absorption |
| Crispness, grain and edge definition of the wash; bare paper | **Inferred uncertainty**: support c |
| The ride: the drop's darkness, pace and sound | **Presentation of measured data**: the observed spectrum played back in order; it adds nothing to it |
| The revision wave after each observation or set-aside | **Presentation** of the change from the old to the new estimate |
| The lifted, hollow stretch of a set-aside ribbon | **Presentation** of the fact that those measured pixels are not used |
| Ink-and-paper material, camera, other sounds, gold and reward effects, ten-drop budget, hints | **Presentation and game mechanics** |

The inferred washes are never presented as observations: they are an estimate from the observed sightlines and change as sightlines are added or set aside.

## 8. Limitations

- **Simulated, not observed.** One PRIYA simulation, one plane, one noise realisation. Real Lyα tomography must also handle continuum fitting, spectral resolution, metal lines, varying noise and irregular sightline positions.
- **Two dimensions only.** A 2D slice, not a 3D reconstruction.
- **Gaussian estimator.** The Lyα forest is not Gaussian: saturated absorption, skewness. The conditional mean is the best linear estimate, not the full posterior.
- **Fixed smoothing.** Structure finer than 1.5 Mpc/h between sightlines is not claimed.
- **Covariance from the same simulation,** estimated from sightlines outside the surveyed region.
- **A simple damped-absorber rule,** a threshold on one spectrum, as described in §5.
- **No comparison with the true field in the game.** The player never sees the true field. Developer checks against it were used only to validate the method.

## 9. Parameters

| Quantity | Value |
|---|---|
| Redshift | z = 3 |
| Plane | 480 sightlines × 1389 pixels (0.25 Mpc/h × 10 km/s); 120 Mpc/h periodic |
| Survey | 25 candidate sightlines, 1.5 Mpc/h apart; 10 observations |
| Shown depth | 40–80 Mpc/h (160 bins of 0.25 Mpc/h) |
| Noise | σ_F = 0.02 per 20 km/s |
| Target smoothing | Gaussian σ = 1.5 Mpc/h, both directions |
| Line-of-sight modes | \|q\| ≤ 60 (all 241 for set-aside conditioning) |
| Covariance | 250 sightlines outside the survey; taper 15–30 Mpc/h |
| Damped-absorber rule | core A > 0.93 for ≥ 40 × 20 km/s; extended while A > 0.3 |
| Support shown as unresolved | c < 0.30 |
