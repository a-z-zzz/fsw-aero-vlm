import os
import sys
import pandas as pd
import numpy as np
from scipy import stats
from scipy.interpolate import interp1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt



try:
    from x29_vlm_validation import (
        AR, E_REF, k_v, build_x29, build_conventional,
        calc_oswald_from_polar, get_spanwise_lift_distribution,
        _extract_vlm_distribution, _interp_cdi_at_cl
    )
except ImportError:
    print("Could not import from x29_vlm_validation. Ensure you run this from the root workspace directory.")
    sys.exit(1)

OUTPUT_DIR = "figures_x29/highres_figures"
os.makedirs(OUTPUT_DIR, exist_ok=True)
LEGEND_DIR = "figures_x29/legends"
os.makedirs(LEGEND_DIR, exist_ok=True)

def extract_and_save_legend(ax, base_filename, info_text=None):
    # Ensure no legend is drawn on the plot
    legend = ax.get_legend()
    if legend is not None:
        legend.remove()

    handles, labels = ax.get_legend_handles_labels()
    
    txt_filename = os.path.join(LEGEND_DIR, f"{base_filename}_legend.txt")
    with open(txt_filename, 'w') as f:
        f.write("Legend entries for this figure:\n")
        f.write("-" * 40 + "\n")
        for label in labels:
            f.write(f"- {label}\n")
            
        if info_text:
            f.write("\nInfo Box / Annotations:\n")
            f.write("-" * 40 + "\n")
            f.write(info_text + "\n")
    return txt_filename

def _plot_polar_validation(df: pd.DataFrame, e_inv: float, e_corr: float, base_name: str):
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
    
    ax.grid(True, alpha=0.3)

    pct_inv = abs(e_inv - E_REF) / E_REF * 100
    pct_corr = abs(e_corr - E_REF) / E_REF * 100
    note = (f'Inviscid Error: {pct_inv:.1f}%\n'
            f'Corrected Error: {pct_corr:.1f}%\n'
            f'Viscous Factor: {k_v:.0%}')

    # Legend handling
    extract_and_save_legend(ax, base_name, info_text=note)

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f"{base_name}.png"), dpi=600, bbox_inches='tight')
    plt.close()

def _plot_lift_curve(df: pd.DataFrame, base_name: str):
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
    
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='k', linewidth=0.5)
    ax.axvline(x=0, color='k', linewidth=0.5)

    note = f'$C_L$$_α$ = {slope:.4f} /deg\n       = {slope * 57.3:.3f} /rad'

    extract_and_save_legend(ax, base_name, info_text=note)

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f"{base_name}.png"), dpi=600, bbox_inches='tight')
    plt.close()

def _plot_spanwise_theoretical(fsw_airplane, asw_airplane, alpha_deg: float, base_name: str):
    fig, ax = plt.subplots(figsize=(12, 8))

    y_fsw, cl_fsw = get_spanwise_lift_distribution(fsw_airplane, alpha_deg)

    ax.plot(y_fsw, cl_fsw, 'b-o', linewidth=2.5, markersize=6,
            label='FSW (X-29, Λ = -29.3°)')

    y_asw = np.linspace(0, 1, 25)
    cl_asw = np.cos(np.pi / 2 * y_asw ** 0.8)
    ax.plot(y_asw, cl_asw, 'r--s', linewidth=2.5, markersize=6,
            label='Conventional ASW (Λ = +29.3°)')

    cl_elliptic = np.sqrt(1 - y_fsw ** 2)
    ax.plot(y_fsw, cl_elliptic, 'k:', linewidth=2, label='Elliptic (Ideal)')

    ax.set_xlabel('Spanwise Station (y / (b/2))', fontsize=14)
    ax.set_ylabel('Normalized Local Lift Coefficient ($C_l / C_{l,max}$)', fontsize=14)
    ax.set_title(f'Spanwise Lift Distribution at α = {alpha_deg}°', fontsize=16)
    
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.2)

    note = ("FSW: Inboard shift reduces tip loading\n"
            "ASW: Outboard loading increases tip stall risk")

    extract_and_save_legend(ax, base_name, info_text=note)

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f"{base_name}.png"), dpi=600, bbox_inches='tight')
    plt.close()

def _plot_spanwise_vlm(fsw_airplane, asw_airplane, alpha_deg: float, base_name: str):
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
    
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 1.05)
    ax.set_ylim(0, 1.15)

    note = ("FSW: Lift centered inboard (Lower root bending stress)\n"
            "ASW: Lift centered outboard (Higher root bending stress)")

    extract_and_save_legend(ax, base_name, info_text=note)

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f"{base_name}.png"), dpi=600, bbox_inches='tight')
    plt.close()

def _plot_ldi_comparison(df_fsw: pd.DataFrame, df_conv: pd.DataFrame, base_name: str,
                          e_corr_fsw: float = None, e_corr_conv: float = None):
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
    
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 16)

    # Use live-computed Oswald efficiencies, not hardcoded constants
    fsw_e_str = f'{e_corr_fsw:.3f}' if e_corr_fsw is not None else 'N/A'
    conv_e_str = f'{e_corr_conv:.3f}' if e_corr_conv is not None else 'N/A'
    note = (f'FSW L/Dᵢ Advantage: +{avg_adv:.1f}%\n'
            f'(Average over α = 4-10° cruise regime)\n\n'
            f'FSW Oswald e = {fsw_e_str}\n'
            f'ASW Oswald e = {conv_e_str}')

    extract_and_save_legend(ax, base_name, info_text=note)

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f"{base_name}.png"), dpi=600, bbox_inches='tight')
    plt.close()

def _plot_total_ld(df_fsw: pd.DataFrame, df_conv: pd.DataFrame, base_name: str, cd0: float = 0.022):
    fig, ax = plt.subplots(figsize=(12, 9))

    df_fp = df_fsw[df_fsw['CL'] > 0.05].copy()
    df_cp = df_conv[df_conv['CL'] > 0.05].copy()

    df_fp['CD_total'] = cd0 + df_fp['CDi_corrected']
    df_fp['L_D_total'] = df_fp['CL'] / df_fp['CD_total']

    df_cp['CD_total'] = cd0 + df_cp['CDi_corrected']
    df_cp['L_D_total'] = df_cp['CL'] / df_cp['CD_total']

    ax.plot(df_fp['CL'], df_fp['L_D_total'], 'b-o', linewidth=2.5, markersize=8,
            label='FSW (X-29, Λ = -29.3°)')
    ax.plot(df_cp['CL'], df_cp['L_D_total'], 'r--s', linewidth=2.5, markersize=8,
            label='Conventional ASW (Λ = +29.3°)')

    ax.set_xlabel('Lift Coefficient ($C_L$)', fontsize=14)
    ax.set_ylabel('Total Aerodynamic Efficiency ($L/D$)', fontsize=14)
    ax.set_title(f'Total L/D vs Lift Coefficient (Estimated $C_{{D0}} = {cd0}$)', fontsize=16)
    
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, max(df_fp['CL'].max(), df_cp['CL'].max()) * 1.05)

    note = (f'Note: This is an estimation.\n'
            f'Zero-lift drag $C_{{D0}}={cd0}$ applied uniformly.\n'
            f'Plotted against $C_L$ to isolate aerodynamic efficiency.')

    extract_and_save_legend(ax, base_name, info_text=note)

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f"{base_name}.png"), dpi=600, bbox_inches='tight')
    plt.close()


def _plot_polar_fsw_asw(df_fsw: pd.DataFrame, df_conv: pd.DataFrame, base_name: str,
                         e_corr_fsw: float = None, e_corr_conv: float = None):
    fig, ax = plt.subplots(figsize=(12, 9))

    df_fp = df_fsw[df_fsw['CL'] > 0.05].copy()
    df_cp = df_conv[df_conv['CL'] > 0.05].copy()

    # Use live-computed Oswald efficiencies in legend labels
    fsw_e_str = f'{e_corr_fsw:.3f}' if e_corr_fsw is not None else '?'
    conv_e_str = f'{e_corr_conv:.3f}' if e_corr_conv is not None else '?'

    ax.plot(df_fp['CL_squared'], df_fp['CDi_corrected'],
            color='blue', linestyle='-', linewidth=3,
            marker='o', markersize=8, markerfacecolor='blue', markeredgecolor='white',
            label=f'FSW (X-29, Λ = -29.3°) — e = {fsw_e_str}', zorder=3)
    ax.plot(df_cp['CL_squared'], df_cp['CDi_corrected'],
            color='red', linestyle='--', linewidth=3,
            marker='s', markersize=8, markerfacecolor='red', markeredgecolor='white',
            label=f'ASW (Λ = +29.3°) — e = {conv_e_str}', zorder=2)

    cdi_fsw_05 = _interp_cdi_at_cl(df_fsw, 0.5)
    cdi_asw_05 = _interp_cdi_at_cl(df_conv, 0.5)

    note = ""
    if not np.isnan(cdi_fsw_05) and not np.isnan(cdi_asw_05):
        cl2_05 = 0.25
        ax.axvline(x=cl2_05, color='gray', linestyle=':', linewidth=1.5, alpha=0.7)
        ax.scatter([cl2_05], [cdi_fsw_05], color='blue', s=150,
                   marker='*', zorder=5, edgecolors='white', linewidths=2)
        ax.scatter([cl2_05], [cdi_asw_05], color='red', s=150,
                   marker='*', zorder=5, edgecolors='white', linewidths=2)

        drag_diff = ((cdi_asw_05 - cdi_fsw_05) / cdi_asw_05) * 100
        note = (f'At exactly CL = 0.5:\n'
                f'  FSW CDi = {cdi_fsw_05:.5f}\n'
                f'  ASW CDi = {cdi_asw_05:.5f}\n'
                f'  FSW drag advantage: {drag_diff:.2f}%')

    ax.set_xlabel('$C_L^2$', fontsize=14)
    ax.set_ylabel('$C_{Di}$ (Induced Drag Coefficient)', fontsize=14)
    ax.set_title('Induced Drag Polar: FSW vs Conventional Swept Wing', fontsize=16)
    
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, df_fp['CL_squared'].max() * 1.05)
    ax.set_ylim(0, max(df_cp['CDi_corrected'].max(), df_fp['CDi_corrected'].max()) * 1.25)

    extract_and_save_legend(ax, base_name, info_text=note)

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f"{base_name}.png"), dpi=600, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    print("Loading CSV data...")
    csv_dir = "data"
    
    df_fsw = pd.read_csv(os.path.join(csv_dir, 'x29_vlm_drag_polar_mar22.csv'))
    df_conv = pd.read_csv(os.path.join(csv_dir, 'reg_asw_vlm_drag_polar_mar22.csv'))

    # Compute Oswald efficiency for both configurations from the live CSV data
    e_inv_fsw, e_corr_fsw, _ = calc_oswald_from_polar(df_fsw)
    e_inv_conv, e_corr_conv, _ = calc_oswald_from_polar(df_conv)

    print(f"  FSW  e_inv={e_inv_fsw:.4f}  e_corr={e_corr_fsw:.4f}")
    print(f"  ASW  e_inv={e_inv_conv:.4f}  e_corr={e_corr_conv:.4f}")

    print("Generating pure CSV-based plots at 600 DPI (no legends, no info boxes)...")
    _plot_polar_validation(df_fsw, e_inv_fsw, e_corr_fsw, "fig1_drag_polar_validation_HighRes")
    _plot_lift_curve(df_fsw, "fig2_lift_curve_HighRes")
    _plot_total_ld(df_fsw, df_conv, "fig6_total_LD_ratio_HighRes")
    _plot_ldi_comparison(df_fsw, df_conv, "fig4_fsw_vs_conventional_HighRes",
                          e_corr_fsw=e_corr_fsw, e_corr_conv=e_corr_conv)
    _plot_polar_fsw_asw(df_fsw, df_conv, "fig5_drag_polar_comparison_HighRes",
                         e_corr_fsw=e_corr_fsw, e_corr_conv=e_corr_conv)

    # ── Spanwise distribution plots commented out ──────────────────────────
    # These require a fresh single-point VLM solve and are not part of the
    # journal submission figures (fig 1, 5, 6). Uncomment to regenerate locally.
    #
    # print("Building geometries and running single-point VLM for spanwise plots...")
    # x29_ap = build_x29()
    # conv_ap = build_conventional()
    # _plot_spanwise_theoretical(x29_ap, conv_ap, 10.0, "fig3_spanwise_lift_distribution_HighRes")
    # _plot_spanwise_vlm(x29_ap, conv_ap, 10.0, "fig3_real_vlm_spanwise_HighRes")

    print(f"All done! Check the '{OUTPUT_DIR}' folder for your images and legend text files.")
