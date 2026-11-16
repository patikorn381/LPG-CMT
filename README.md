# LPG-CMT: transmission spectra of long-period fiber gratings and cladding-mode selection

Open MATLAB code (with independent Python checks) for the paper

> S. Sutharak and P. Anchuen, "Simulation of power transmission spectra in long period grating optical fiber and cladding-mode selection with a turning-point-independent figure of merit."

The code computes the effective indices of the LP01 core mode and the LP0m cladding modes of a step-index fiber, the coupling coefficients, and the transmission spectrum of a uniform LPG with coupled-mode theory (scalar, weak guidance). It also computes the SRI response, the figure of merit FOM = |S|/FWHM, and the validation checks reported in the paper.

## Requirements
- MATLAB R2021a or later (no toolboxes required)
- Optional: Python 3 with numpy, scipy, pandas and matplotlib, for the independent checks and Figure 8

## How to run
Open MATLAB in this folder, then run:

| Command | Produces | Paper |
|---|---|---|
| `run_all` | Tables 1–4, Figures 1–3, convergence, material dispersion and tolerance checks; log in `results/summary.txt`. Also writes supplementary figures not shown in the paper (mode profiles, coupling coefficients, shift versus SRI) | Methods, Results 1–2 |
| `added_analyses` | DC self-coupling, band depth, FWHM and FOM versus SRI | Results 1.6, Table 5 |
| `upper_branch_nu15` | Second (upper-branch) band of ν = 15 versus SRI | Results 1.2 |
| `fom_map` | 2,100-design FOM map (`results/fom_map.csv`) and Table 6 (`results/table7_design.csv`) | Methods 2.4, Results 2.3 |
| `bench_libish` | Comparison with the measured LPGs of Libish et al. (`results/bench_libish.csv`) | Table 3 |

Optional Python checks (run in this folder):

| Command | Produces | Paper |
|---|---|---|
| `python verify_independent.py` | Nine numerical verification checks | Results 1.4 |
| `python vector_check.py` | Comparison with the exact three-layer vector model (HE/EH, l = 1) | Tables 3–4 |
| `python make_fig8.py` | `results/fig8_fom_map.png` from `results/fom_map.csv` | Figure 4 |

All parameters are set in `lpg_params.m` (Table 1 of the paper).

## Main files
| File | Purpose |
|---|---|
| `lpg_params.m` | Fiber and grating parameters |
| `neff_core.m`, `disp2_core.m` | LP01 effective index, two-layer model (Eq. 1) |
| `neff_clad.m`, `disp3_clad.m` | LP0m cladding-mode effective indices, three-layer model |
| `bisect.m` | Bisection root finder (tolerance 1e-12) |
| `coupling.m` | Coupling coefficient from the overlap integral (Eq. 3) |
| `lpg_spectrum.m`, `band_fwhm.m` | Transmission spectrum, band position, depth and FWHM (Eqs. 2, 4) |
| `sri_sweep.m` | Resonant wavelength versus SRI |
| `xhalf.m` | Half-width parameter x_h(κL) used in Eqs. (6)–(7) |
| `confinement.m` | Core power fractions for the DC self-coupling term |
| `index_model.m`, `sellmeier_silica.m` | Constant index (default) or silica Sellmeier dispersion |
| `lpg_main.m` | Main calculation and Figures 1–7 |
| `validation_specs.m`, `compare_singh.m` | NA, cut-off, mode-field diameter, and literature comparisons |
| `dispersion_check.m`, `convergence_check.m`, `tolerance_analysis.m` | Material dispersion, grid convergence and tolerance checks |

## Output
- `results/` holds the CSV and text outputs used in the paper, plus `fig8_fom_map.png` (Figure 4 of the paper). The other figures are written to the same folder when `run_all` is run.

## Model assumptions
1. Scalar (LP) approximation. The (m−1)-th cladding mode LP0m corresponds to the odd-order mode ν = 2m − 3 of Erdogan's vector model (ν = 15 is LP09).
2. Uniform grating, cosine index modulation Δn uniform over the core. The DC self-coupling term is omitted in the main model and assessed separately (`added_analyses.m`).
3. Each cladding mode is treated as a two-mode coupled problem, and the total spectrum is the product of the individual responses.
4. Material dispersion is omitted in the main model and assessed separately (`dispersion_check.m`).

## License
Released for research and teaching use. Please cite the paper above if you use this code.
