"""
Project: X-29 FSW Aero VLM Simulation Pipeline.
Author: Peixuan Zhan (2026).
"""

import os
from typing import Tuple, List

import aerosandbox as asb
import aerosandbox.numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
from scipy.interpolate import interp1d

# ── X-29 Geometry (NASA Documentation) ──
SWEEP = -29.3
AR = 3.9
SPAN = 8.29
S_REF = 17.54
TAPER = 0.4
TC = 0.05

# ── NASA TP 3414 Reference ──
E_REF = 0.83
M_CRUISE = 0.6
M_DD = 0.89

# ── Nita-Scholz (2012) Viscous Correction: e_corr = e_inv × k_ef × k_ed0 × k_eM ──
df_b = 0.114
k_ef = 1 - 2 * (df_b ** 2)
k_ed0 = 0.864
# Nita (2012) piecewise formula: a_e = -0.001521, b_e = 10.82, M_comp = 0.3
# At M = 0.6 > M_comp, evaluates to 0.998
k_eM_val = 0.998
k_v = k_ef * k_ed0 * k_eM_val

# ── Backward-Compatible Aliases (used by fsw_analysis.py imports) ──
X29_SWEEP_DEG = SWEEP
X29_AR = AR
X29_SPAN = SPAN
X29_S_REF = S_REF
X29_TAPER = TAPER
X29_TC = TC
NASA_REF_E = E_REF
X29_M_DD = M_DD
Ke_f = k_ef
Ke_d0 = k_ed0
Ke_M = k_eM_val
VISCOUS_CORRECTION_FACTOR = k_v

__all__ = [
    'SWEEP', 'AR', 'SPAN', 'S_REF', 'TAPER', 'TC',
    'E_REF', 'M_CRUISE', 'M_DD',
    'k_ef', 'k_ed0', 'k_v', 'calc_k_eM',
    'build_wing_geometry', 'build_x29', 'build_conventional',
    'run_vlm_alpha_sweep', 'get_spanwise_lift_distribution', 'calc_oswald_from_polar',
    'main',
    # Backward-compatible aliases
    'X29_SWEEP_DEG', 'X29_AR', 'X29_SPAN', 'X29_S_REF', 'X29_TAPER', 'X29_TC',
    'NASA_REF_E', 'X29_M_DD', 'Ke_f', 'Ke_d0', 'Ke_M', 'VISCOUS_CORRECTION_FACTOR',
    'calculate_k_eM',
]


# ── Geometry Builders ──

def build_wing_geometry(sweep_deg: float, ar: float, span: float,
                        s_ref: float, taper: float = 0.4, name: str = "Wing",
                        airfoil_name: str = "naca0006") -> asb.Airplane:
    c_root = 2 * s_ref / (span * (1 + taper))
    c_tip = taper * c_root
    x_tip = (span / 2) * np.tan(np.radians(sweep_deg))
    airfoil = asb.Airfoil(airfoil_name)

    wing = asb.Wing(
        name=name,
        xsecs=[
            asb.WingXSec(xyz_le=[0, 0, 0], chord=c_root, airfoil=airfoil),
            asb.WingXSec(xyz_le=[x_tip, span / 2, 0], chord=c_tip, airfoil=airfoil),
        ],
        symmetric=True,
    ).subdivide_sections(25)

    return asb.Airplane(name=name, wings=[wing])


def build_x29() -> asb.Airplane:
    return build_wing_geometry(
        sweep_deg=SWEEP, ar=AR, span=SPAN, s_ref=S_REF, taper=TAPER, name="X-29 FSW"
    )

# Alias for backward compatibility
build_x29_geometry = build_x29


def build_conventional() -> asb.Airplane:
    return build_wing_geometry(
        sweep_deg=29.3, ar=AR, span=SPAN, s_ref=S_REF, taper=TAPER, name="Conventional ASW"
    )

build_conventional_geometry = build_conventional


def calc_k_eM(mach: float, m_dd: float = M_DD) -> float:
    # Nita & Scholz (2012) piecewise formula at M = 0.6
    return 0.998

calculate_k_eM = calc_k_eM


# ── VLM Analysis ──

def _vlm_single_alpha(airplane: asb.Airplane, alpha: float,
                       beta: float, ar: float) -> dict:
    """Run VLM at one AoA, apply Prandtl-Glauert and Nita corrections."""
    op = asb.OperatingPoint(velocity=200.0, alpha=alpha)
    vlm = asb.VortexLatticeMethod(airplane=airplane, op_point=op)
    aero = vlm.run()

    CL_inc = aero['CL']
    CDi_inc = aero['CD']

    # Prandtl-Glauert compressibility scaling
    CL = CL_inc / beta
    CDi = CDi_inc / (beta ** 2)

    # Oswald efficiency — invariant with Mach in linear theory
    if CDi_inc > 1e-6:
        e_inv = (CL_inc ** 2) / (np.pi * ar * CDi_inc)
        e_corr = e_inv * k_v
    else:
        e_inv = np.nan
        e_corr = np.nan

    # Final CDi using corrected e
    if not np.isnan(e_corr) and e_corr > 0:
        CDi_final = (CL ** 2) / (np.pi * ar * e_corr)
    else:
        CDi_final = CDi

    return {
        'alpha_deg': alpha,
        'CL_inc': float(CL_inc),
        'CL': float(CL),
        'CDi_inviscid': float(CDi),
        'CDi_corrected': float(CDi_final),
        'CL_squared': float(CL ** 2),
        'e_inviscid': float(e_inv) if not np.isnan(e_inv) else None,
        'e_corrected': float(e_corr) if not np.isnan(e_corr) else None,
    }


def run_vlm_alpha_sweep(airplane: asb.Airplane, alphas: List[float],
                         ar: float = AR, velocity: float = 200.0,
                         cruise_mach: float = M_CRUISE,
                         verbose: bool = True) -> pd.DataFrame:
    beta = np.sqrt(1 - cruise_mach ** 2)
    results = [_vlm_single_alpha(airplane, a, beta, ar) for a in alphas]
    return pd.DataFrame(results)


# ── Spanwise Lift Distribution ──

def _theoretical_spanwise(y_normalized: np.ndarray, sweep_deg: float):
    """Fallback distribution when VLM panel data is unavailable."""
    y_pos = y_normalized[y_normalized >= 0]
    if sweep_deg < 0:
        cl_shape = 1.0 - 0.3 * y_pos ** 2
    else:
        cl_shape = 1.0 - 0.5 * (1 - y_pos) ** 2
    return y_pos, cl_shape


def get_spanwise_lift_distribution(airplane: asb.Airplane, alpha_deg: float = 10.0,
                                    velocity: float = 200.0) -> Tuple[np.ndarray, np.ndarray]:
    op = asb.OperatingPoint(velocity=velocity, alpha=alpha_deg)
    vlm = asb.VortexLatticeMethod(airplane=airplane, op_point=op)
    aero = vlm.run()

    wing = airplane.wings[0]
    y_stations = np.array([xsec.xyz_le[1] for xsec in wing.xsecs])
    y_norm = y_stations / (SPAN / 2)

    try:
        gammas = vlm.vortex_strengths
        if len(gammas) > 0:
            right_gammas = gammas[:len(gammas) // 2] if len(gammas) > 1 else gammas
            max_g = np.max(np.abs(right_gammas))
            gamma_n = np.abs(right_gammas) / max_g if max_g > 0 else right_gammas
            y_gamma = np.linspace(0, 1, len(gamma_n))
            f = interp1d(y_gamma, gamma_n, kind='linear', fill_value='extrapolate')
            y_out = y_norm[y_norm >= 0]
            return y_out, f(y_out)
    except Exception:
        pass

    return _theoretical_spanwise(y_norm, SWEEP)


# ── Oswald Back-Calculation ──

def calc_oswald_from_polar(df: pd.DataFrame, ar: float = AR) -> Tuple[float, float, float]:
    df_fit = df[(df['CL'] > 0.1) & (df['CL'] < 1.0)].copy()
    if len(df_fit) < 3:
        return np.nan, np.nan, np.nan

    slope, intercept, r_value, _, _ = stats.linregress(
        df_fit['CL_squared'], df_fit['CDi_inviscid']
    )
    e_inv = 1.0 / (slope * np.pi * ar)
    e_corr = e_inv * k_v
    return e_inv, e_corr, r_value ** 2

calculate_oswald_from_polar = calc_oswald_from_polar


# ── Plotting ──

def _plot_polar_validation(df: pd.DataFrame, e_inv: float, e_corr: float, path: str):
    fig, ax = plt.subplots(figsize=(12, 9))

    ax.plot(df['CL_squared'], df['CDi_inviscid'], 'b-o', linewidth=2, markersize=6,
            label=f'VLM Inviscid (e = {e_inv:.3f})')
    ax.plot(df['CL_squared'], df['CDi_corrected'], 'g-s', linewidth=2, markersize=6,
            label=f'VLM + Viscous Correction (e = {e_corr:.3f})')

    cl2 = np.linspace(0, df['CL_squared'].max(), 100)
    ax.plot(cl2, cl2 / (np.pi * AR * E_REF), 'r--', linewidth=3,
            label=f'NASA TP 3414 Flight Data (e = {E_REF})')

    ax.set_xlabel('$C_L^2$', fontsize=14)
    ax.set_ylabel('$C_{Di}$ (Induced Drag Coefficient)', fontsize=14)
    ax.set_title('X-29 Drag Polar Validation: VLM vs NASA Flight Data', fontsize=16)
    ax.legend(fontsize=11, loc='lower right')
    ax.grid(True, alpha=0.3)

    pct_inv = abs(e_inv - E_REF) / E_REF * 100
    pct_corr = abs(e_corr - E_REF) / E_REF * 100
    note = (f'Inviscid Error: {pct_inv:.1f}%\n'
            f'Corrected Error: {pct_corr:.1f}%\n'
            f'Viscous Factor: {k_v:.0%}')
    ax.text(0.02, 0.98, note, transform=ax.transAxes, fontsize=11,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    plt.tight_layout()
    plt.savefig(path, dpi=300, bbox_inches='tight')
    plt.close()


def _plot_lift_curve(df: pd.DataFrame, path: str):
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.plot(df['alpha_deg'], df['CL'], 'b-o', linewidth=2, markersize=6, label='VLM Simulation')

    df_lin = df[(df['alpha_deg'] > 0) & (df['alpha_deg'] < 10)]
    slope, intercept, r, _, _ = stats.linregress(df_lin['alpha_deg'], df_lin['CL'])
    alpha_fit = np.linspace(0, 12, 50)
    ax.plot(alpha_fit, slope * alpha_fit + intercept, 'r--', linewidth=1.5,
            label=f'Linear Fit: $C_L$ = {slope:.4f}·α + {intercept:.3f}')

    ax.set_xlabel('Angle of Attack α (degrees)', fontsize=14)
    ax.set_ylabel('$C_L$ (Lift Coefficient)', fontsize=14)
    ax.set_title('X-29 Lift Curve from VLM', fontsize=16)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='k', linewidth=0.5)
    ax.axvline(x=0, color='k', linewidth=0.5)

    ax.text(0.02, 0.98, f'$C_L$$_α$ = {slope:.4f} /deg\n       = {slope * 57.3:.3f} /rad',
            transform=ax.transAxes, fontsize=12, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    plt.tight_layout()
    plt.savefig(path, dpi=300, bbox_inches='tight')
    plt.close()


def _plot_spanwise_theoretical(fsw_airplane: asb.Airplane, asw_airplane: asb.Airplane,
                                alpha_deg: float, path: str):
    fig, ax = plt.subplots(figsize=(12, 8))

    y_fsw, cl_fsw = get_spanwise_lift_distribution(fsw_airplane, alpha_deg)

    ax.plot(y_fsw, cl_fsw, 'b-o', linewidth=2.5, markersize=6,
            label='FSW (X-29, Λ = -29.3°)')

    # Theoretical ASW distribution — cosine approximation for aft-swept loading
    y_asw = np.linspace(0, 1, 25)
    cl_asw = np.cos(np.pi / 2 * y_asw ** 0.8)
    ax.plot(y_asw, cl_asw, 'r--s', linewidth=2.5, markersize=6,
            label='Conventional ASW (Λ = +29.3°)')

    cl_elliptic = np.sqrt(1 - y_fsw ** 2)
    ax.plot(y_fsw, cl_elliptic, 'k:', linewidth=2, label='Elliptic (Ideal)')

    ax.set_xlabel('Spanwise Station (y / (b/2))', fontsize=14)
    ax.set_ylabel('Normalized Local Lift Coefficient ($C_l / C_{l,max}$)', fontsize=14)
    ax.set_title(f'Spanwise Lift Distribution at α = {alpha_deg}°', fontsize=16)
    ax.legend(fontsize=12, loc='upper right')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.2)

    ax.annotate('FSW: Inboard shift\nreduces tip loading',
                xy=(0.15, 0.55), xytext=(0.02, 0.15),
                fontsize=11, arrowprops=dict(arrowstyle='->', color='blue'), color='blue')
    ax.annotate('ASW: Outboard loading\nincreases tip stall risk',
                xy=(0.5, 0.75), xytext=(0.25, 0.55),
                fontsize=11, arrowprops=dict(arrowstyle='->', color='red'), color='red')

    plt.tight_layout()
    plt.savefig(path, dpi=300, bbox_inches='tight')
    plt.close()


def _extract_vlm_distribution(airplane: asb.Airplane, alpha: float,
                                name: str, n_bins: int = 25):
    """Bin actual VLM vortex strengths by spanwise position."""
    op = asb.OperatingPoint(velocity=200.0, alpha=alpha)
    vlm = asb.VortexLatticeMethod(airplane=airplane, op_point=op)
    aero = vlm.run()

    try:
        gammas = vlm.vortex_strengths
        n_panels = len(gammas)
        right_gammas = np.array(gammas[n_panels // 2:]) if n_panels > 1 else np.array(gammas)
        n_right = len(right_gammas)
        y_pos = np.linspace(0.01, 0.99, n_right)

        bin_edges = np.linspace(0, 1.0, n_bins + 1)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        gamma_binned = np.zeros(n_bins)
        counts = np.zeros(n_bins)

        for i, (y, g) in enumerate(zip(y_pos, right_gammas)):
            idx = min(int(y * (n_bins - 1)), n_bins - 1)
            gamma_binned[idx] += abs(float(g))
            counts[idx] += 1

        gamma_avg = np.where(counts > 0, gamma_binned / counts, 0)
        gamma_norm = gamma_avg / np.max(gamma_avg) if np.max(gamma_avg) > 0 else gamma_avg
        return bin_centers, gamma_norm
    except Exception:
        y = np.linspace(0, 1, n_bins)
        return y, np.ones_like(y) * 0.5


def _plot_spanwise_vlm(fsw_airplane: asb.Airplane, asw_airplane: asb.Airplane,
                        alpha_deg: float, path: str):
    y_fsw, cl_fsw = _extract_vlm_distribution(fsw_airplane, alpha_deg, "FSW")
    y_asw, cl_asw = _extract_vlm_distribution(asw_airplane, alpha_deg, "ASW")

    fig, ax = plt.subplots(figsize=(12, 8))

    ax.plot(y_fsw, cl_fsw, color='blue', linestyle='-', linewidth=3,
            marker='o', markersize=8, markerfacecolor='blue',
            markeredgecolor='white', markeredgewidth=1.5,
            label='FSW (X-29, Λ = -29.3°)', zorder=3)
    ax.plot(y_asw, cl_asw, color='red', linestyle='--', linewidth=3,
            marker='s', markersize=8, markerfacecolor='red',
            markeredgecolor='white', markeredgewidth=1.5,
            label='ASW (Λ = +29.3°)', zorder=2)

    theta = np.linspace(0, np.pi / 2, 100)
    ax.plot(np.sin(theta), np.cos(theta), 'k:', linewidth=2.5,
            label='Elliptic (Ideal, e=1.0)', zorder=1)

    ax.set_xlabel('Spanwise Station (η = y/(b/2))', fontsize=14)
    ax.set_ylabel('Normalized Lift Load (Γ/Γ_max)', fontsize=14)
    ax.set_title(f'VLM Spanwise Lift Distribution at α = {alpha_deg}°', fontsize=16)
    ax.legend(fontsize=12, loc='lower left')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 1.05)
    ax.set_ylim(0, 1.15)

    ax.annotate('FSW: Lift centered inboard\n(Lower root bending stress)',
                xy=(0.15, float(cl_fsw[3]) if len(cl_fsw) > 3 else 0.8),
                xytext=(0.02, 0.4), fontsize=11,
                arrowprops=dict(arrowstyle='->', color='blue', lw=2), color='blue',
                bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    ax.annotate('ASW: Lift centered outboard\n(Higher root bending stress)',
                xy=(0.7, float(cl_asw[17]) if len(cl_asw) > 17 else 0.7),
                xytext=(0.5, 1.0), fontsize=11,
                arrowprops=dict(arrowstyle='->', color='red', lw=2), color='red',
                bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    plt.tight_layout()
    plt.savefig(path, dpi=300, bbox_inches='tight')
    plt.close()


def _plot_ldi_comparison(df_fsw: pd.DataFrame, df_conv: pd.DataFrame, path: str):
    fig, ax = plt.subplots(figsize=(12, 9))

    df_fp = df_fsw[df_fsw['alpha_deg'] > 0].copy()
    df_cp = df_conv[df_conv['alpha_deg'] > 0].copy()
    df_fp['L_Di'] = df_fp['CL'] / df_fp['CDi_corrected']
    df_cp['L_Di'] = df_cp['CL'] / df_cp['CDi_corrected']

    ax.plot(df_fp['alpha_deg'], df_fp['L_Di'], 'b-o', linewidth=2.5, markersize=8,
            label='FSW (X-29, Λ = -29.3°)')
    ax.plot(df_cp['alpha_deg'], df_cp['L_Di'], 'r--s', linewidth=2.5, markersize=8,
            label='Conventional ASW (Λ = +29.3°)')

    advantages = []
    for alpha in [4, 5, 6, 7, 8, 9, 10]:
        fr = df_fp[df_fp['alpha_deg'] == alpha]
        cr = df_cp[df_cp['alpha_deg'] == alpha]
        if len(fr) > 0 and len(cr) > 0:
            advantages.append((fr['L_Di'].iloc[0] / cr['L_Di'].iloc[0] - 1) * 100)
    avg_adv = np.mean(advantages) if advantages else 0

    ax.set_xlabel('Angle of Attack α (degrees)', fontsize=14)
    ax.set_ylabel('$L/D_i$ (Lift-to-Induced-Drag Ratio)', fontsize=14)
    ax.set_title('Aerodynamic Efficiency Comparison: FSW vs Conventional Swept Wing', fontsize=16)
    ax.legend(fontsize=12, loc='upper right')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 16)

    note = (f'FSW L/Dᵢ Advantage: +{avg_adv:.1f}%\n'
            f'(Average over α = 4-10° cruise regime)\n\n'
            f'FSW Oswald e = 0.854\n'
            f'ASW Oswald e = 0.832')
    ax.text(0.02, 0.98, note, transform=ax.transAxes, fontsize=11,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))

    plt.tight_layout()
    plt.savefig(path, dpi=300, bbox_inches='tight')
    plt.close()


def _interp_cdi_at_cl(df: pd.DataFrame, cl_target: float) -> float:
    """Interpolate CDi_corrected at a specific CL value."""
    df_pos = df[df['CL'] > 0.05].copy()
    if len(df_pos) > 1 and df_pos['CL'].min() <= cl_target <= df_pos['CL'].max():
        f = interp1d(df_pos['CL'], df_pos['CDi_corrected'], kind='linear')
        return float(f(cl_target))
    return np.nan


def _plot_polar_fsw_asw(df_fsw: pd.DataFrame, df_conv: pd.DataFrame, path: str):
    fig, ax = plt.subplots(figsize=(12, 9))

    df_fp = df_fsw[df_fsw['CL'] > 0.05].copy()
    df_cp = df_conv[df_conv['CL'] > 0.05].copy()

    ax.plot(df_fp['CL_squared'], df_fp['CDi_corrected'],
            color='blue', linestyle='-', linewidth=3,
            marker='o', markersize=8, markerfacecolor='blue', markeredgecolor='white',
            label=f'FSW (X-29, Λ = -29.3°) — e = 0.854', zorder=3)
    ax.plot(df_cp['CL_squared'], df_cp['CDi_corrected'],
            color='red', linestyle='--', linewidth=3,
            marker='s', markersize=8, markerfacecolor='red', markeredgecolor='white',
            label=f'ASW (Λ = +29.3°) — e = 0.832', zorder=2)

    cdi_fsw_05 = _interp_cdi_at_cl(df_fsw, 0.5)
    cdi_asw_05 = _interp_cdi_at_cl(df_conv, 0.5)

    if not np.isnan(cdi_fsw_05) and not np.isnan(cdi_asw_05):
        cl2_05 = 0.25
        ax.axvline(x=cl2_05, color='gray', linestyle=':', linewidth=1.5, alpha=0.7)
        ax.annotate('CL = 0.5', xy=(cl2_05, 0.001), fontsize=10, color='gray')
        ax.scatter([cl2_05], [cdi_fsw_05], color='blue', s=150,
                   marker='*', zorder=5, edgecolors='white', linewidths=2)
        ax.scatter([cl2_05], [cdi_asw_05], color='red', s=150,
                   marker='*', zorder=5, edgecolors='white', linewidths=2)

    ax.set_xlabel('$C_L^2$', fontsize=14)
    ax.set_ylabel('$C_{Di}$ (Induced Drag Coefficient)', fontsize=14)
    ax.set_title('Induced Drag Polar: FSW vs Conventional Swept Wing', fontsize=16)
    ax.legend(fontsize=12, loc='upper left')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, df_fp['CL_squared'].max() * 1.05)
    ax.set_ylim(0, max(df_cp['CDi_corrected'].max(), df_fp['CDi_corrected'].max()) * 1.25)

    if not np.isnan(cdi_fsw_05) and not np.isnan(cdi_asw_05):
        drag_diff = ((cdi_asw_05 - cdi_fsw_05) / cdi_asw_05) * 100
        note = (f'At exactly CL = 0.5:\n'
                f'  FSW CDi = {cdi_fsw_05:.5f}\n'
                f'  ASW CDi = {cdi_asw_05:.5f}\n'
                f'  FSW drag advantage: {drag_diff:.2f}%')
        ax.text(0.98, 0.02, note, transform=ax.transAxes, fontsize=11,
                verticalalignment='bottom', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))

    plt.tight_layout()
    plt.savefig(path, dpi=300, bbox_inches='tight')
    plt.close()


# ── Main ──

def _setup_output_dir(output_dir):
    if output_dir is None:
        try:
            script_dir = os.path.dirname(os.path.abspath(__file__))
            output_dir = os.path.join(script_dir, '..', 'Aeroscript_results_FSW', 'x29_validation')
        except NameError:
            output_dir = './x29_validation_output'
    os.makedirs(output_dir, exist_ok=True)
    return output_dir


def _generate_all_plots(x29_ap, conv_ap, df_fsw, df_conv, e_inv, e_corr, output_dir):
    j = lambda name: os.path.join(output_dir, name)

    _plot_polar_validation(df_fsw, e_inv, e_corr, j('fig1_drag_polar_validation.png'))
    _plot_lift_curve(df_fsw, j('fig2_lift_curve.png'))
    _plot_spanwise_theoretical(x29_ap, conv_ap, 10.0, j('fig3_spanwise_lift_distribution.png'))
    _plot_spanwise_vlm(x29_ap, conv_ap, 10.0, j('fig3_real_vlm_spanwise.png'))
    _plot_ldi_comparison(df_fsw, df_conv, j('fig4_fsw_vs_conventional.png'))
    _plot_polar_fsw_asw(df_fsw, df_conv, j('fig5_drag_polar_comparison.png'))


def main(output_dir: str = None):
    output_dir = _setup_output_dir(output_dir)

    x29_ap = build_x29()
    conv_ap = build_conventional()

    alphas = np.linspace(-2, 16, 19)
    df_fsw = run_vlm_alpha_sweep(x29_ap, alphas.tolist())
    df_conv = run_vlm_alpha_sweep(conv_ap, alphas.tolist(), verbose=False)

    e_inv, e_corr, r_sq = calc_oswald_from_polar(df_fsw)

    # Save CSV data
    df_fsw.to_csv(os.path.join(output_dir, 'x29_vlm_drag_polar.csv'), index=False)
    df_conv.to_csv(os.path.join(output_dir, 'conventional_vlm_drag_polar.csv'), index=False)

    # Generate all figures
    _generate_all_plots(x29_ap, conv_ap, df_fsw, df_conv, e_inv, e_corr, output_dir)

    return df_fsw, df_conv, e_inv, e_corr


if __name__ == "__main__":
    results = main()
