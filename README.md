# Ink Forest

**An interactive Lyα tomography game.**

You have ten sightlines into a universe you cannot see. Ride each Lyα forest spectrum, reconstruct the structure between them, and decide where to look next. Some sightlines reveal the cosmic web; others can distort the map.

**Play:** https://jibancat.github.io/ink-forest/ (any recent desktop, tablet or phone browser with WebGL2; sound on).

## What you do

The game is one loop: **observe** a sightline, **infer** the structure between your sightlines, **revise** the map and your next choice.

A short rules card opens the game; press **Begin**.

1. **Choose a sightline.** Each light above the sheet is the line of sight to a distant quasar. Click it, or on a touch screen tap it twice.
2. **Ride it.** A drop of ink falls down the sightline. It slows and darkens where intergalactic gas absorbs the quasar's light. What you see and hear is that sightline's Lyα forest spectrum, and only the part already passed. On the way the drop collects coins: one for each 1 Mpc/h-wide strip of the sheet this line will newly resolve, so lines far from the ones you know pay most. Now and then one is a diamond, worth 5.
3. **Watch the map change.** When the ride ends, the sheet between your sightlines is inferred again from everything you have observed. More coins appear wherever the map changed.
4. **Decide where to look next.** You have ten drops, and the map can never be complete.
5. **Set a line aside if it looks wrong.** Press and hold a landed line to set aside its damped stretch; hold it again to restore it.

At the end the map holds still. Then a short debrief shows your gold and how much of the sheet you resolved, and separates what you measured from what you inferred.

## The science in brief

**A Lyα forest sightline.** Light from a distant quasar crosses billions of light-years of intergalactic hydrogen. Neutral hydrogen absorbs at the Lyman-α wavelength. Because the universe expands, gas at each distance absorbs at a different observed wavelength. The spectrum therefore becomes a one-dimensional map of the gas along that line of sight: a "forest" of absorption lines.

**Measured, inferred, presentation.**

| On screen | Category | What it is |
|---|---|---|
| Ink ribbons | **Measured** | The spectrum of each sightline you observed, drawn from that line's own (simulated, noisy) data at full resolution. |
| Washes between the ribbons | **Inferred** | A statistical estimate of the large-scale gas distribution, computed only from the sightlines you observed. Darker means more absorbing gas than average (denser regions); paler means less (emptier regions). |
| How settled the ink looks; bare paper | **Inferred uncertainty** | Crisp, dry ink is well constrained; grainy, wet ink is uncertain; bare paper is not constrained by any sightline. |
| The ride, the ink-and-paper look, camera, sounds, the ten-drop budget | **Presentation and game mechanics** | The ride plays the measured spectrum back in order (its pace, darkness and sound follow it) but adds no information to it. The rest is presentation. |
| Coins and diamonds | **Game rewards** | Ride coins count the new ground a line will resolve, from sightline positions only; more coins mark where the map changed; diamonds are luck. They never depend on how much a line absorbs and never mark dense gas ([SCIENCE.md](SCIENCE.md) §7). |

**The map is an estimate, never the answer.** The game never shows the simulation's true structure, and the washes are not observed: they are inferred from your sightlines and change as you add more.

**Why a damped absorber can distort the map.** Occasionally a sightline passes through a single dense cloud of gas, a damped Lyα absorber. That cloud absorbs all the light over a long stretch of the spectrum. Read as ordinary forest, it looks like a huge wall of gas, and the reconstruction spreads that wall onto the neighbouring sightlines. Setting the stretch aside removes those pixels from the inference. The map then withdraws to what the neighbours actually support, and that region becomes less certain rather than being filled in. The game never tells you which line is which.

## What is simulated, and what is not claimed

- **The spectra are simulated,** from the PRIYA cosmological hydrodynamic simulations, not observed with a telescope. Each spectrum carries added noise at a fixed level. No instrument, continuum or contaminating metal lines are modelled.
- **The map is one two-dimensional slice,** 37.5 × 40 Mpc/h at redshift z = 3. It is neither a three-dimensional reconstruction nor a reconstruction of the cosmic web as a whole.
- **The reconstruction is deliberately simple:** a Gaussian conditional mean at a fixed 1.5 Mpc/h smoothing. It shows how Lyα tomography works; it is not a survey-grade pipeline.

[SCIENCE.md](SCIENCE.md) explains exactly what each element of the picture means.

## Run it yourself

Static files only: no build step, no backend, no accounts, no tracking or analytics, and no external fonts or libraries.

```
python3 -m http.server 8000
```

Then open http://localhost:8000/. Opening `index.html` directly from disk also works in most browsers.

## Data and licence

- **Code and documentation:** MIT licence, © 2026 Ming-Feng Ho ([LICENSE](LICENSE)).
- **Data:** the PRIYA-derived runtime data in `assets/ink-forest-data.js` are not covered by the MIT licence. They are redistributed with permission for this project. See [DATA.md](DATA.md) for contents, provenance and terms.

## Citation

See [CITATION.cff](CITATION.cff). If you use the data, please also cite the PRIYA simulations:

> S. Bird, M. Fernandez, M.-F. Ho, M. Qezlou, R. Monadi, Y. Ni, N. Chen, R. Croft, T. Di Matteo, *PRIYA: A New Suite of Lyman-alpha Forest Simulations for Cosmology*, JCAP 10 (2023) 037, doi:10.1088/1475-7516/2023/10/037, arXiv:2306.05471.
