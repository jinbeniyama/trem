#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot corner map with TI of rocks, TI of regolith, Hapke angle, and alpha.
Note: Assume that bond albedo A is fixed.

"""
from argparse import ArgumentParser as ap
import pandas as pd
import numpy as np

from trem.emittance.common_emittance import calc_confidence_chi2


if __name__ == "__main__":
    parser = ap(description="Plot a corner plot.")
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

    import corner
    param_cols = ['Htheta', 'TIrego', 'TIrock', 'alpha']

    
    # Choose reliable once
    df = df[df["chi2"] < chi2_min + chi2_3sigma]

    # Remove specific columns to avoid a following error
    # > ValueError: It looks like the parameter(s) in column(s) 1 have no dynamic range. Please provide a `range` argument.
    for p in param_cols:
        N_p = len(set(df[p]))
        if N_p == 1:
            param_cols.remove(p)
            print(f"Not plot {p}")

    data_array = df[param_cols].values
    fig = corner.corner(
        data_array, labels=param_cols, 
        #show_titles=True, 
        label_kwargs={"fontsize": 10}, 
        #title_kwargs={"fontsize": 12}, 
        smoonth=1,
        plot_datapoints=True,      # サンプル点も表示
        plot_density=True,         # 密度塗りつぶし
        plot_contours=True,
        bins=50
        )
    fig.savefig(args.out, dpi=300, bbox_inches='tight')
