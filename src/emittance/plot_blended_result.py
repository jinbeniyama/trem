#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot 3-d color-map with (x, y, z) = (TI of rocks, TI of regolith, Hapke angle)
Note: Assume that bond albedo A is fixed.

Example
-------
# Plot 3-d chi2 map.
> plot_blended_result.py result.txt
"""
from argparse import ArgumentParser as ap
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from trem.emittance.common_emittance import calc_confidence_chi2


if __name__ == "__main__":
    parser = ap(description="Plot 3-d chi2 map after blending fluxes.")
    parser.add_argument(
        "res", type=str, 
        help="Results of blending.")
    parser.add_argument(
        "--TI_thresh", type=float, default=1500,
        help="Threshold of TI of regolith and TI of rocks")
    parser.add_argument(
        "--dof", type=int, default=1,
        help="Defree of freedom")
    parser.add_argument(
        "--paper", type=str, default="P14",
        help="P14 or V17, type of uncertainty")
    parser.add_argument(
        "--reduce", action="store_true", default=False,
        help="Reduced chi square")
    parser.add_argument(
        "--out", type=str, default="blend_chi2_map.jpg",
        help="Output file")
    args = parser.parse_args()
   
    
    # Read results
    df = pd.read_csv(args.res, sep=" ")
    N_all = len(df)
    title = f"TI_th = {args.TI_thresh}"

    dof = args.dof
    # Calculate reduced chi2 
    if args.reduce:
        print(f"  Calculate reduced chi square with dof = {dof}")
        df["chi2"] = df["chi2"]/dof

    # Extract minimum chi2 and its index
    idx_min = df["chi2"].idxmin()
    chi2_min = df.loc[idx_min, "chi2"]
    TIrego_min = df.loc[idx_min, "TIrego"]
    TIrock_min = df.loc[idx_min, "TIrock"]
    Htheta_min = df.loc[idx_min, "Htheta"]
    alpha_min = df.loc[idx_min, "alpha"]

    # Add 1-sigma, 3-sigma
    chi2_1sigma = calc_confidence_chi2(args.paper, chi2_min, dof, 1, args.reduce)
    chi2_3sigma = calc_confidence_chi2(args.paper, chi2_min, dof, 3, args.reduce)

    # Uncertainties of TIs 
    chi2_arr = np.array(df["chi2"])
    #   TI of regolith with chi2 < chi2_min + chi2_3sigma
    TIrego_arr = np.array(df["TIrego"])
    TIrego_arr_sig = TIrego_arr[chi2_arr < chi2_min + chi2_3sigma]
    TIrego3sigl, TIrego3sigu = np.min(TIrego_arr_sig), np.max(TIrego_arr_sig)
    #   TI of rock with chi2 < chi2_min + chi2_3sigma
    TIrock_arr = np.array(df["TIrock"])
    TIrock_arr_sig = TIrock_arr[chi2_arr < chi2_min + chi2_3sigma]
    TIrock3sigl, TIrock3sigu = np.min(TIrock_arr_sig), np.max(TIrock_arr_sig)
    #   Htheta with chi2 < chi2_min + chi2_3sigma
    Htheta_arr = np.array(df["Htheta"])
    Htheta_arr_sig = Htheta_arr[chi2_arr < chi2_min + chi2_3sigma]
    Htheta3sigl, Htheta3sigu = np.min(Htheta_arr_sig), np.max(Htheta_arr_sig)
    #   alpha with chi2 < chi2_min + chi2_3sigma
    alpha_arr = np.array(df["alpha"])
    alpha_arr_sig = alpha_arr[chi2_arr < chi2_min + chi2_3sigma]
    alpha3sigl, alpha3sigu = np.min(alpha_arr_sig), np.max(alpha_arr_sig)

    text = (
        f"TIrego = ${TIrego_min}_" + "{" + f"-{TIrego_min-TIrego3sigl}" + "}^" 
        "{" + f"+{TIrego3sigu-TIrego_min}" + "}" + f"$ (N={len(TIrego_arr_sig)})\n"
        f"TIrock = ${TIrock_min}_" + "{" + f"-{TIrock_min-TIrock3sigl}" + "}^" 
        "{" + f"+{TIrock3sigu-TIrock_min}" + "}" + f"$ (N={len(TIrock_arr_sig)})\n"
        f"Htheta = ${Htheta_min}_" + "{" + f"-{Htheta_min-Htheta3sigl}" + "}^" 
        "{" + f"+{Htheta3sigu-Htheta_min}" + "}" + f"$ (N={len(Htheta_arr_sig)})\n"
        f"alpha = ${alpha_min}_" + "{" + f"-{alpha_min-alpha3sigl:.2f}" + "}^" 
        "{" + f"+{alpha3sigu-alpha_min:.2f}" + "}" + f"$ (N={len(alpha_arr_sig)})\n"
        )

    key_x, key_y, key_z, key_c = "TIrego", "TIrock", "Htheta", "chi2"
    xlabel = "TI of regolith [tiu]"
    ylabel = "TI of rocks [tiu]"
    zlabel = "Hapke angle [deg]"

    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_axes([0.05, 0.10, 0.7, 0.7], projection = '3d')
    cax = fig.add_axes([0.85, 0.10, 0.02, 0.7])
    ax.grid(which="both", color="gray",linewidth=0.2)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)

    # Add range of TIrego and TIrock
    ax.text(0.1, 0.1, 1.0, text, size=12, transform=ax.transAxes)

    # Histogram of alpha
    ax_a = fig.add_axes([0.82, 0.88, 0.10, 0.06])
    ax_a.set_xlabel(r"$\alpha$", fontsize=12)
    ax_a.set_ylabel("Fraction", fontsize=12)

    
    x = df[key_x]
    y = df[key_y]
    z = df[key_z]
    c = df[key_c]
    # Add minimum value

    label_min = (
            f"Minimum chi2 {chi2_min:.2f} at (TIrego, TIrock, Htheta, alpha) = ({TIrego_min:.1f}, {TIrock_min:.1f}, {Htheta_min:.1f}, {alpha_min:.2f})")
    ax.scatter(
        df.loc[idx_min, key_x], df.loc[idx_min, key_y], df.loc[idx_min, key_z], 
        color="red", marker="x", s=100, label=label_min)
    pl = ax.scatter(
        x, y, z, c=c, cmap='viridis', label=f"N={N_all}")
    cbar = fig.colorbar(pl, cax=cax)
    cbar.set_label("$\chi^2$")

    # Histogram of alpha
    for a in sorted(list(set(df["alpha"]))):
        df_temp = df[df["alpha"] == a]
        print(f"a = {a:.2f} N = {len(df_temp)}")
    ax_a.hist(df["alpha"])
    ax_a.set_yscale("log")
    ymin, ymax = ax_a.get_ylim()
    ax_a.set_ylim([1, ymax])

    ax.legend(bbox_to_anchor=(1.05, 1), borderaxespad=0)
    fig.savefig(args.out)
