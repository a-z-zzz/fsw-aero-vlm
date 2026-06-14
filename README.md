# fsw-aero-vlm
A non-linear Vortex Lattice Method (VLM) Python simulation pipeline for evaluating forward-swept wing (FSW) aerodynamics. Used in conjunction with a physics-based NASA X-29 validation and trade study tooling utilizing AeroSandbox.

# X-29 FSW VLM Simulation Pipeline
[![DOI](https://img.shields.io/badge/DOI-pending-lightgrey.svg)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## Abstract
The Vortex Lattice Method (VLM), a lightweight computational model, offers faster alternatives for aeronautical modeling than Computational Fluid Dynamics (CFD). However, their accuracy for unconventional designs, like Forward Swept Wings (FSW), is not well-developed. We assess the capability of using a Python-based VLM solver to simulate static aerodynamic behavior of FSW planforms by using the Grumman X-29 as a test case. We simulate the main control surfaces of the Grumman X-29, apply Prandtl-Glauert and viscous corrections, and compare our simulated drag polars and Oswald efficiency factor to published NASA flight data. A geometrically mirrored Aft-Swept Wing (ASW) equivalent is also simulated to isolate theoretical sweep effects. We find that our simulated Oswald efficiency factor (0.743) is 10.5% lower than the NASA flight data (~0.83). However, the solver correctly finds that the forward sweep does not directly alter L/D or induced drag compared to an equivalent ASW under inviscid conditions. VLM can, in summary, approximate general aerodynamic trends, but requires design-specific calibrations to produce more quantitatively accurate findings.


---

## Methdology

The Vortex Lattice Method approximates the lifting surface as a discrete array of panels, each with a horseshoe vortex whose strength Γ is determined by enforcing a no-penetration boundary condition at the panel control point. The velocity induced by each vortex is computed via the Biot-Savart law, allowing the full spanwise downwash distribution to be resolved across all panels. By the Kutta-Joukowski theorem, the local aerodynamic force per unit span is L' = ρV∞Γ, which, when decomposed into freestream-aligned components, returns lift and induced drag. The latter arises from the rear tilt of the local vector by the induced angle αᵢ = w/V∞, giving Di' ≈ L' × αᵢ for the small angles typical of cruising flight conditions. Total forces are recovered by spanwise integration across all panel contributions, and the Oswald efficiency factor is back-solved from the resulting CL and CDi values using the Prandtl drag polar. Because VLM assumes inviscid, incompressible potential flow, a viscous correction is applied in post-processing following the method of Nita and Scholz (2018).

---

## Key Results

| Metric | Value |
|--------|-------|
| Validation dataset | NASA TP 3414 (Grumman X-29, Mach 0.6) |
| VLM Inviscid span efficiency | $e_{inv} = 0.872$ |
| VLM Corrected span efficiency | $e_{corr} = 0.733$ |
| Inviscid validation error | 5.1% against NASA reference ($e = 0.83$) |
| Corrected validation error | 11.7% against NASA reference |
| Viscous correction factor $k_v$ | 84% (Nita–Scholz method) |
| FSW Oswald efficiency | $e = 0.854$ |
| ASW Oswald efficiency | $e = 0.832$ |
| Induced drag difference (FSW vs ASW at $C_L = 0.5$) | –0.02% (marginal) |
| Total Aerodynamic Efficiency ($L/D$) | FSW sustains higher maximum $L/D$ across cruise $C_L$ (see Fig. 6) |

> **Note on validation error:** The 11.7% corrected span efficiency error reflects the Nita–Scholz viscous correction being applied to the wing-only inviscid VLM model. The full configuration (with canard and strakes) is expected to reduce this discrepancy. The inviscid model agrees to within 5.1%.

---

## Installation & Dependencies

Python 3.10 or higher is required. All core dependencies can be installed via pip:
```bash
pip install aerosandbox numpy scipy matplotlib pandas
```

### Optional: OpenVSP & SU2

For high-fidelity mesh generation and RANS verification, ensure `OpenVSP 3.x` and `SU2 v8.0` are available in your system path. These are not required to run the VLM pipeline but are used for higher-fidelity cross-validation.

- [OpenVSP Download](http://openvsp.org/)
- [SU2 Download](https://su2code.github.io/)

---

## Usage & Execution

### 1. X-29 VLM Validation (Full Configuration)

Verifies solver accuracy against NASA TP 3414 flight data utilizing the complete X-29 geometry, incorporating the closely-coupled canard and leading-edge strakes to capture complex upwash/downwash interactions:

[bash]
python3 x29_vlm_validation.py
[/bash]

Expected outputs:
- Total aerodynamic efficiency ($L/D$) comparison vs lift coefficient
- Coupled spanwise lift distributions accounting for canard downwash
- Drag polar overlay against NASA reference data
- Console report of $e_{inviscid}$, $k_v$, $e_{corrected}$, and validation error

### 2. Boeing 737-800 Trade Study

Runs the aerodynamic trade study comparing FSW and ASW configurations of identical planform geometry:
```bash
python3 Aeroscript_FSW/fsw_analysis.py --num_sims 100
```

Expected outputs:
- Induced drag polar comparison at cruise ($C_L = 0.5$)
- Oswald efficiency summary for FSW and ASW configurations
- Aerodynamic efficiency comparison plot ($L/D_i$ vs $ lpha$)

---

### 4. Generate High-Res and Journal Plots

The visual results for publication can be programmatically reproduced by running the export scripts after the main simulation:
```bash
python3 export_journal_plots.py
python3 export_high_res_plots.py
```
This ensures reproducibility without manual tweaking, generating clean, monochromatic figures for print and high-fidelity colored figures for digital use.

---

## Repository Structure
```
fsw-aero-vlm/
├── x29_vlm_validation.py              # Main X-29 solver validation against NASA TP 3414
├── export_journal_plots.py            # Reproduces monochromatic journal figures
├── export_high_res_plots.py           # Reproduces high-resolution color figures
├── data/
│   ├── x29_vlm_drag_polar_mar22.csv      # FSW Drag Polar raw numbers
│   └── reg_asw_vlm_drag_polar_mar22.csv  # ASW Drag Polar raw numbers
├── figures_x29/                     
│   ├── journal_figures/               # Output for journal publication
│   └── legends/                       # Raw text legends
└── README.md
```

---

## Citation

If you use this code in your research, please cite the repository:
```bibtex
@misc{zhan2026fswvlm,
  author       = {Zhan, Andy},
  title        = {fsw-aero-vlm},
  year         = {2026},
  doi          = {10.5281/zenodo.XXXXXXX},
  note         = {DOI pending. Repository available at https://github.com/a-z-zzz/fsw-aero-vlm}
}
```
## License

This project is licensed under the MIT License. See [`LICENSE`](LICENSE) for details.
