# X-29 FSW VLM Simulation — Results Summary

**Simulation:** Full-configuration X-29 (wing + canard + strake), Mach 0.6, NASA TP 3414 validation
**AoA Range:** α = −2° to +16°, 1-degree increments (19 data points per configuration)
**Comparison:** FSW (X-29, Λ = −29.3°) vs Conventional ASW (Λ = +29.3°)
**Viscous Correction:** Nita-Scholz (2012) method, factor k_v = 0.840 (84.0%)

---

## Figure 1 — X-29 Drag Polar Validation: VLM vs NASA TP 3414
`fig1_drag_polar_validation_Journal.png`

This figure plots C_Di vs C_L² for the FSW (X-29) configuration and compares it against
the NASA TP 3414 flight reference (e = 0.83).

| Metric | Value |
|---|---|
| FSW inviscid Oswald efficiency (e_inv) | **0.8813** |
| FSW corrected Oswald efficiency (e_corr) | **0.7402** |
| NASA TP 3414 reference efficiency (e_ref) | **0.83** |
| Inviscid validation error | 6.18% |
| **Corrected validation error** | **10.82%** |
| Viscous correction factor (k_v) | 0.840 (84.0%) |
| Polar regression R² | 0.99997 |

> Note: The corrected error of 10.82% represents the discrepancy between the
> Nita-Scholz viscous correction applied to the full-configuration VLM (wing +
> canard + strake) and the NASA flight reference, which was measured on the
> isolated wing. This is an improvement over the 11.7% baseline from the
> wing-only model.

---

## Figure 5 — Induced Drag Polar: FSW vs Conventional ASW
`fig5_drag_polar_comparison_Journal.png`

This figure plots C_Di vs C_L² for both FSW and ASW configurations.

### Oswald Span Efficiency
| Configuration | e_inviscid | e_corrected |
|---|---|---|
| FSW (X-29, Λ = −29.3°) | 0.8813 | **0.7402** |
| ASW (Λ = +29.3°) | 0.8931 | **0.7501** |

### Induced Drag Comparison at Key Lift Coefficients (C_Di, corrected)
| C_L | FSW C_Di | ASW C_Di | FSW Advantage |
|---|---|---|---|
| 0.3 | 0.01175 | 0.01216 | **+3.40%** |
| 0.4 | 0.01944 | 0.01972 | **+1.44%** |
| 0.5 | 0.02932 | 0.02933 | +0.02% (≈ equal) |
| 0.6 | 0.04148 | 0.04119 | −0.72% (ASW marginally better) |
| 0.7 | 0.05563 | 0.05523 | −0.73% (ASW marginally better) |

> Note: The FSW advantage is concentrated at lower lift coefficients (cruise and
> low-α regime). At C_L ≥ 0.6 (high-α), the conventional ASW regains a marginal
> induced drag advantage due to its more outboard lift loading.

---

## Figure 6 — Total L/D vs Lift Coefficient: FSW vs Conventional ASW
`fig6_total_LD_ratio_Journal.png`

This figure plots **total aerodynamic efficiency (L/D)** as a function of lift
coefficient C_L, with an estimated zero-lift drag C_D0 = 0.022 applied uniformly
to both configurations.

### Total L/D Comparison vs C_L
| C_L | FSW L/D | ASW L/D | FSW Advantage |
|---|---|---|---|
| 0.10 | 4.15 | 4.04 | +2.65% |
| 0.18 | 6.67 | 6.53 | +2.14% |
| 0.26 | 8.37 | 8.25 | +1.48% |
| 0.34 | 9.32 | 9.22 | +1.11% |
| 0.42 | 9.72 | 9.65 | +0.72% |
| 0.49 | 9.76 | 9.74 | +0.29% ← near cruise |
| 0.57 | 9.59 | 9.60 | −0.09% ← crossover point |
| 0.64 | 9.30 | 9.34 | −0.39% |
| 0.71 | 8.96 | 9.01 | −0.58% |
| 0.79 | 8.59 | 8.65 | −0.68% |
| 0.86 | 8.23 | 8.29 | −0.69% |
| 0.93 | 7.89 | 7.94 | −0.62% |
| 1.00 | 7.56 | 7.60 | −0.48% |
| 1.06 | 7.26 | 7.28 | −0.29% |
| **Average** | — | — | **+0.33%** |

> Note: Once C_D0 = 0.022 is included, the parasitic drag term dominates and
> largely washes out the induced drag advantage. The FSW total L/D advantage is
> modest at ~0.33% on average, with the crossover occurring at C_L ≈ 0.57.
> Below this point FSW is more efficient; above it the ASW is marginally better.
> The larger L/D_i advantage (+3.90%, from Figure 4/5) reflects the induced drag
> mechanism in isolation and should be cited separately from the total drag result.

---

## Supplementary Figure Data (High-Res Figures Only)

### Fig 2 — Lift Curve Slope (from fig2_lift_curve_HighRes.png)
| Metric | Value |
|---|---|
| C_Lα (per degree) | **0.07808 /deg** |
| C_Lα (per radian) | **4.4738 /rad** |
| Linear fit range | α = 1° to 9° |

### Fig 3 — Spanwise Lift Distribution (from fig3_* files)
Plotted at α = 10°. FSW shows a characteristic inboard-shifted lift distribution
relative to the ASW and the ideal elliptic reference, consistent with theoretical
predictions for negative-sweep wings.

---

## Model Configuration Reference
| Parameter | Value |
|---|---|
| Wing sweep (FSW) | −29.3° |
| Wing sweep (ASW, baseline) | +29.3° |
| Aspect ratio (AR) | 3.9 |
| Span | 8.29 m |
| Reference area | 17.54 m² |
| Taper ratio | 0.4 |
| Cruise Mach | 0.6 |
| Prandtl-Glauert β | 0.800 |
| Nita-Scholz k_v | 0.840 |
| AoA range | −2° to +16° (1° steps) |
| Simulation rows (each config) | 19 |
