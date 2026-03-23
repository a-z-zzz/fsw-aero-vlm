# fsw-aero-vlm
A non-linear Vortex Lattice Method (VLM) Python simulation pipeline for evaluating forward-swept wing (FSW) aerodynamics. Used in conjunction with a physics-based NASA X-29 validation and trade study tooling utilizing AeroSandbox.

# Boeing 737-800 FSW Simulation Pipeline
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A validated Vortex Lattice Method (VLM) solver and simulation suite for auditing the aerodynamic feasibility of Forward-Swept Wing (FSW) configurations on commercial transport airframes.

---

## Abstract
This technical note investigates the aerodynamic feasibility of forward-swept wing (FSW) configurations for a Boeing 737-800 class aircraft using a validated Vortex Lattice Method (VLM). The solver was calibrated against Grumman X-29 flight test data, yielding a corrected span efficiency of $e = 0.851$ within 2.56% of the NASA TP 3414 reference value. At cruise conditions ($C_L = 0.5$), the FSW configuration produces an Oswald efficiency of $e = 0.854$ compared to $e = 0.832$ for a geometrically identical aft-swept baseline. 

---

## Key Results
* **Validation:** VLM calibrated against NASA TP 3414 X-29 flight data at Mach 0.6.
* **Accuracy:** Resolved Oswald efficiency ($e_{corrected} = 0.851$) with a marginal 2.56% error.
* **Efficiency Delta:** FSW configurations demonstrated a **3.14% reduction in induced drag** (~1.26% total drag reduction) over geometrically mirrored aft-swept baselines.

---

## Installation & Dependencies

### Pre-requisites
To run this pipeline, you must have **Python 3.10+** installed. Computational modules depend on the following libraries:

```bash
pip install numpy scipy matplotlib pandas
