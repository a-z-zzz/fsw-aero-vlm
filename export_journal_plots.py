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
        AR, E_REF, calc_oswald_from_polar, _interp_cdi_at_cl
    )
except ImportError:
    print("Could not import from x29_vlm_validation.")
    sys.exit(1)

OUTPUT_DIR = "Journal_Export"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Set matplotlib to use a monochromatic style for legends
plt.rcParams['legend.frameon'] = True
plt.rcParams['legend.facecolor'] = 'white'
plt.rcParams['legend.edgecolor'] = 'black'


def _plot_polar_validation(df: pd.DataFrame, e_inv: float, e_corr: float, base_name: str):
    fig, ax = plt.subplots(figsize=(10, 8))

    ax.plot(df['CL_squared'], df['CDi_inviscid'], 'b-o', linewidth=2, markersize=6,
            label='VLM Inviscid')
    ax.plot(df['CL_squared'], df['CDi_corrected'], 'g-s', linewidth=2, markersize=6, fillstyle='none',
            label='VLM Viscous Corrected')

    cl2 = np.linspace(0, df['CL_squared'].max(), 100)
    ax.plot(cl2, cl2 / (np.pi * AR * E_REF), 'r--', linewidth=2,
            label='NASA Flight Data')

    ax.set_xlabel('$C_L^2$', fontsize=14)
    ax.set_ylabel('Induced Drag Coefficient $C_{Di}$', fontsize=14)
    ax.set_title('Drag Polar Validation', fontsize=16)
    


    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=12, loc='lower right', framealpha=1.0)

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f"{base_name}.png"), dpi=600, bbox_inches='tight')
    plt.close()

def _plot_polar_fsw_asw(df_fsw: pd.DataFrame, df_conv: pd.DataFrame, base_name: str):
    fig, ax = plt.subplots(figsize=(10, 8))

    df_fp = df_fsw[df_fsw['CL'] > 0.05].copy()
    df_cp = df_conv[df_conv['CL'] > 0.05].copy()

    ax.plot(df_fp['CL_squared'], df_fp['CDi_corrected'],
            color='blue', linestyle='-', linewidth=2,
            marker='o', markersize=6,
            label='FSW Configuration', zorder=3)
    ax.plot(df_cp['CL_squared'], df_cp['CDi_corrected'],
            color='red', linestyle='--', linewidth=2,
            marker='s', markersize=6, fillstyle='none',
            label='ASW Configuration', zorder=2)

    ax.set_xlabel('$C_L^2$', fontsize=14)
    ax.set_ylabel('Induced Drag Coefficient $C_{Di}$', fontsize=14)
    ax.set_title('Induced Drag Polar', fontsize=16)
    


    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, df_fp['CL_squared'].max() * 1.05)
    ax.set_ylim(0, max(df_cp['CDi_corrected'].max(), df_fp['CDi_corrected'].max()) * 1.25)

    ax.legend(fontsize=12, loc='upper left', framealpha=1.0)

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f"{base_name}.png"), dpi=600, bbox_inches='tight')
    plt.close()

def _plot_total_ld(df_fsw: pd.DataFrame, df_conv: pd.DataFrame, base_name: str, cd0: float = 0.022):
    fig, ax = plt.subplots(figsize=(10, 8))

    df_fp = df_fsw[df_fsw['CL'] > 0.05].copy()
    df_cp = df_conv[df_conv['CL'] > 0.05].copy()

    df_fp['CD_total'] = cd0 + df_fp['CDi_corrected']
    df_fp['L_D_total'] = df_fp['CL'] / df_fp['CD_total']

    df_cp['CD_total'] = cd0 + df_cp['CDi_corrected']
    df_cp['L_D_total'] = df_cp['CL'] / df_cp['CD_total']

    ax.plot(df_fp['CL'], df_fp['L_D_total'], 'b-o', linewidth=2, markersize=6,
            label='FSW Configuration')
    ax.plot(df_cp['CL'], df_cp['L_D_total'], 'r--s', linewidth=2, markersize=6, fillstyle='none',
            label='ASW Configuration')

    ax.set_xlabel('Lift Coefficient $C_L$', fontsize=14)
    ax.set_ylabel('Total L/D Ratio', fontsize=14)
    ax.set_title('Total Aerodynamic Efficiency', fontsize=16)
    


    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, max(df_fp['CL'].max(), df_cp['CL'].max()) * 1.05)

    ax.legend(fontsize=12, loc='lower right', framealpha=1.0)

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f"{base_name}.png"), dpi=600, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    print("Loading CSV data...")
    csv_dir = "Aeroscript_results_FSW/x29_validation_final"
    if not os.path.exists(csv_dir):
        csv_dir = "Aeroscript_results_FSW/x29_validation"
        
    df_fsw = pd.read_csv(os.path.join(csv_dir, 'x29_vlm_drag_polar.csv'))
    df_conv = pd.read_csv(os.path.join(csv_dir, 'conventional_vlm_drag_polar.csv'))

    e_inv, e_corr, _ = calc_oswald_from_polar(df_fsw)

    print("Generating Journal Plots (fig 1, 5, 6) in 'Journal_Export'...")
    _plot_polar_validation(df_fsw, e_inv, e_corr, "fig1_drag_polar_validation_Journal")
    _plot_polar_fsw_asw(df_fsw, df_conv, "fig5_drag_polar_comparison_Journal")
    _plot_total_ld(df_fsw, df_conv, "fig6_total_LD_ratio_Journal")

    print(f"All done! Check the '{OUTPUT_DIR}' folder.")
