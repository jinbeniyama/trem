#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plot TIrock vs. TIrego vs. Htheta.
"""
from argparse import ArgumentParser as ap
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd


def plot_TIrego_vs_TIrock(df, TIth, H_list, out="out.jpg", vmin=None, vmax=None, logx=False, logy=False):

    # Make plot
    # Nx = 11 (alpha)
    # Ny = 60 (Htheta)
    nrows, ncols = 11, 60

    TImin, TImax = np.min(df["TIrego"])*0.9, np.max(df["TIrock"])*1.1
    if TImin <= 0:
        TImin = 1e-2

    # chi2 min and max
    vmin = vmin if vmin is not None else np.min(df["chi2"])
    vmax = vmax if vmax is not None else np.max(df["chi2"])

    
    fig = plt.figure(figsize=(80, 25))
    gs = gridspec.GridSpec(nrows, ncols + 1, width_ratios=[1] * ncols + [0.05])
    axes = []


    # i: Index of alpha
    alpha_list = sorted(list(set(df["alpha"])))
    for i in range(nrows):
        alpha_i = alpha_list[i]
        df_a = df[df["alpha"] == alpha_i]
        row_axes = []
        # j: Index of Htheta
        for j in range(ncols):
            H_j = H_list[j]
            df_H_a = df_a[df_a["Htheta"] == H_j]

            ax = fig.add_subplot(gs[i, j]) 
            ax.set_xlim([TImin, TIth*1.1])
            ax.set_ylim([TIth*0.9, TImax])
            if logx:
                ax.set_xscale("log")
            if logy:
                ax.set_yscale("log")
            

            # 軸ラベル表示を減らす・フォント小さく
            if j == 0:
                ax.set_ylabel(r"$\Gamma_{rock}$", fontsize=8)
                ax.text(0.05, 0.85, r"$\alpha$=" + f"{alpha_i:.1f}", fontsize=6, transform=ax.transAxes)
            else:
                ax.set_yticklabels([])

            if i == nrows - 1:
                ax.set_xlabel(r"$\Gamma_{regolith}$", fontsize=8)
            else:
                ax.set_xticklabels([])

            if i == 0:
                ax.set_title(f"{H_j:.1f}", fontsize=6)

            if len(df_H_a) == 0:
                ax.text(0.05, 0.5, "No data", fontsize=6, transform=ax.transAxes)
            else:
                pl = ax.scatter(
                    df_H_a["TIrego"],
                    df_H_a["TIrock"],
                    c=df_H_a["chi2"],
                    s=10,
                    cmap="viridis",
                    alpha=0.4,
                    vmin=vmin,
                    vmax=vmax
                )

            row_axes.append(ax)
        axes.append(row_axes)

    # Colorbar
    cbar_ax = fig.add_subplot(gs[:, -1])
    cb = fig.colorbar(pl, cax=cbar_ax)
    cb.set_label(r"$\chi_\nu^2$", fontsize=12)
    
    plt.tight_layout()
    fig.subplots_adjust(wspace=0.05, hspace=0.05)
    
    if out:
        plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()


if __name__ == "__main__":
    parser = ap(description="Plot TIrock vs. TIrego vs. Htheta.")
    parser.add_argument(
        "res", type=str,
        help="Result of TPM with chi2.")
    parser.add_argument(
        "-x", type=str, default="Htheta", 
        help="x axis of the plot (A, TI, TIrego etc.)")
    parser.add_argument(
        "--dof", type=int, default=1,
        help="Degree of freedom")
    parser.add_argument(
        "--TIth", type=float, default=1,
        help="Threshold of thermal inertia")
    parser.add_argument(
        "--vr", type=float, nargs=2, default=[0, 100],
        help="Value range")
    parser.add_argument(
        "--out", type=str, default="TIrego_TIrock_Htheta.jpg",
        help="Output file")
    args = parser.parse_args()
   
    
    # Read files 
    df = pd.read_csv(args.res, sep=" ")
    dof = args.dof
    
    # reduce
    df["chi2"] = df["chi2"]/dof
    chi2_min = np.min(df["chi2"])
    print(f"reduced chi2 minimum = {chi2_min:.3f}")
    
    # Htheta
    H_list = sorted(list(set(df["Htheta"])))
    vmin, vmax = args.vr
    plot_TIrego_vs_TIrock(df, args.TIth, H_list, args.out, vmin=vmin, vmax=vmax, logy=True)
