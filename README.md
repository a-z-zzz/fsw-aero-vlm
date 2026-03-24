# fsw-aero-vlm
A non-linear Vortex Lattice Method (VLM) Python simulation pipeline for evaluating forward-swept wing (FSW) aerodynamics. Used in conjunction with a physics-based NASA X-29 validation and trade study tooling utilizing AeroSandbox.

# X-29 FSW VLM Simulation Pipeline
[![DOI](https://img.shields.io/badge/DOI-pending-lightgrey.svg)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## Abstract
This technical note investigates the aerodynamic feasibility of simulating forward-swept wing (FSW) configurations using a Vortex Lattice Method (VLM). The solver was calibrated against Grumman X-29 flight-test data, yielding a corrected span efficiency of e = 0.851, within 2.56% of the NASA TP 3414 reference value. At cruise conditions (C_L= 0.5), the FSW configuration produces an Oswald efficiency of e = 0.854 compared to e = 0.832 for a geometrically identical aft-swept baseline, corresponding to a 3.14% reduction in induced drag. Since induced drag accounts for approximately 40% of total cruise drag, this yields a net drag reduction of roughly 1.24%. When combined with estimates of natural laminar flow (NLF) benefits characteristic of forward-swept geometries, this efficiency gain suggests a meaningful improvement in the baseline lift-to-drag ratio for commercial transport applications.

---

## Methdology

The Vortex Lattice Method approximates the lifting surface as a discrete array of panels, each with a horseshoe vortex whose strength Γ is determined by enforcing a no-penetration boundary condition at the panel control point. The velocity induced by each vortex is computed via the Biot-Savart law, allowing the full spanwise downwash distribution to be resolved across all panels. By the Kutta-Joukowski theorem, the local aerodynamic force per unit span is L' = ρV∞Γ, which, when decomposed into freestream-aligned components, returns lift and induced drag. The latter arises from the rear tilt of the local vector by the induced angle αᵢ = w/V∞, giving Di' ≈ L' × αᵢ for the small angles typical of cruising flight conditions. Total forces are recovered by spanwise integration across all panel contributions, and the Oswald efficiency factor is back-solved from the resulting CL and CDi values using the Prandtl drag polar. Because VLM assumes inviscid, incompressible potential flow, a viscous correction is applied in post-processing following the method of Nita and Scholz (2018).

---

## Key Results

| Metric | Value |
|--------|-------|
| Validation dataset | NASA TP 3414 (Grumman X-29, Mach 0.6) |
| Corrected span efficiency | $e_{corrected} = 0.851$ |
| Validation error | 2.56% against NASA reference |
| FSW Oswald efficiency | $e = 0.854$ |
| ASW Oswald efficiency | $e = 0.832$ |
| Induced drag reduction (FSW vs ASW) | 3.14% |
| Total drag reduction | ~1.26% |

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

### 1. NASA X-29 VLM Validation

Verifies solver accuracy against NASA TP 3414 flight data and generates spanwise lift distribution and drag polar comparisons:
```bash
python3 Aeroscript_FSW/x29_vlm_validation.py
```

Expected outputs:
- Spanwise lift distribution comparison (FSW vs ASW vs elliptical ideal)
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
- Aerodynamic efficiency comparison plot (L/D$_i$ vs $\alpha$)

---

## Repository Structure
```
fsw-aero-vlm/
├── x29_vlm_validation.py    # X-29 solver validation against NASA TP 3414
├── data/
│   └── x29_vlm_drag_polar_mar22.csv      # FSW Drag Polar raw numbers
│   └── reg_asw_vlm_drag_polar_mar22.csv  # ASW Drag Polar raw numbers
├── figures_x29/                     # Generated output figures
└── README.md
```

---

## Citation

If you use this code in your research, please cite both the repository and the associated paper:
```bibtex
@misc{zhan2026fswvlm,
  author       = {Zhan, Andy},
  title        = {fsw-aero-vlm: VLM Simulation Pipeline for Forward-Swept Wing Aerodynamic Analysis},
  year         = {2026},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.XXXXXXX},
  note         = {DOI pending. Repository available at https://github.com/a-z-zzz/fsw-aero-vlm}
}
```

---

## License

This project is licensed under the MIT License. See [`LICENSE`](LICENSE) for details.
