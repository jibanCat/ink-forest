# Data and provenance

Ink Forest ships one data file, `assets/ink-forest-data.js` (about 460 KB). It holds simulated Lyα forest spectra derived from the PRIYA simulations. Everything else in this repository is source code and documentation written for Ink Forest. There are no third-party libraries, fonts, images or audio files: all textures and sounds are generated in code.

## Source

**Simulation.** The PRIYA suite of cosmological hydrodynamic simulations for the Lyman-α forest (Bird et al. 2023; full citation below).

**Parent data.** One plane of synthetic Lyα transmitted-flux sightlines from one PRIYA simulation box at redshift z = 3:
- 480 parallel sightlines, 0.25 Mpc/h apart;
- 1389 pixels of 10 km/s each, spanning the 120 Mpc/h periodic box;
- noiseless, at full resolution.

The parent plane is **not** part of this repository. Ink Forest ships only the small derived product described below.

## What the data file contains

Only what the game needs, computed from the parent plane by a fixed script:

| Field | Content | Size |
|---|---|---|
| `obsfull` | Transmitted flux of the **25 candidate sightlines** (every 6th sightline, 1.5 Mpc/h apart), with Gaussian noise added (σ_F = 0.02 per 20 km/s, a fixed seed per sightline), clipped to [−0.05, 1.1], stored as 16-bit integers | 25 × 1389 |
| `fluxfull` | The same noisy flux, clipped to [0, 1], stored as 8-bit integers (used to draw the ribbons) | 25 × 1389 |
| `obs` | The noisy absorption contrast of the 25 sightlines, Fourier-resampled to 0.25 Mpc/h bins, stored as 16-bit integers (the estimator's input) | 25 × 480 |
| `cdd`, `ctd`, `cddfull` | Covariance tables per line-of-sight Fourier mode, as a function of transverse separation, stored as 32-bit floats (the estimator's prior) | 151 lags × 61 or 241 modes |
| `meta` | Scalars: mean absorption, prior variance, noise variance, grid sizes, the sightline indices | — |

**The 25 noisy sightlines** are what the player can observe. They are the only spectra in the file.

**The covariance tables** are the two-point function of the absorption contrast, measured from 250 sightlines of the parent plane outside the survey region:
- symmetrised;
- tapered beyond 15–30 Mpc/h of transverse separation;
- transformed to line-of-sight Fourier modes.

They are summary statistics, not spectra.

## What is not included

- **The noiseless parent plane** (the 480 sightlines above). The raw plane file, `f5.bin`, is not included in this repository.
- **Any hidden true field.** Neither the noiseless absorption field between the sightlines nor any smoothed version of it is included, so nothing in the file reveals the structure the game asks you to infer.
- **Any absorber catalogue.** Damped absorbers are identified at run time from the noisy spectra by the rule in [SCIENCE.md](SCIENCE.md).
- **Any other spectra.** Sightlines other than the 25 candidates enter only through the covariance statistics above.
- **Any other simulation product:** no other planes or boxes, snapshots, particle data or parameter sets.

## Terms

The derived runtime data included with Ink Forest (`assets/ink-forest-data.js`) are redistributed for this public educational project with permission. They are not covered by the software's MIT licence. This permission extends only to this file as part of Ink Forest. It is not a licence for the PRIYA simulations or any of their other data products. See the PRIYA citation below for the underlying simulation work.

## Citation

If you use these data, please cite the PRIYA simulations:

> S. Bird, M. Fernandez, M.-F. Ho, M. Qezlou, R. Monadi, Y. Ni, N. Chen, R. Croft, T. Di Matteo, *PRIYA: A New Suite of Lyman-alpha Forest Simulations for Cosmology*, Journal of Cosmology and Astroparticle Physics 10 (2023) 037. doi:10.1088/1475-7516/2023/10/037. arXiv:2306.05471.
