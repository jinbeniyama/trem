#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot corner map with TI of rocks, TI of regolith, Hapke angle, and alpha.
Note: Assume that bond albedo A is fixed.

"""
from argparse import ArgumentParser as ap
import pandas as pd
import numpy as np
import corner

from trem.emittance.common_emittance import calc_confidence_chi2


def extract_npercent(samples, percent=68):
    """Obtian median and uncertainties.

    Parameters
    ----------
    samples : array-like
        samples
    percent : float
        percentage of interest

    Return
    ------
    med, val_l, val_u : float
        median, lower and upper bounds
    """
    samples = np.asarray(samples)
    med = np.median(samples)

    alpha = (100 - percent) / 2
    lower, upper = np.percentile(samples, [alpha, 100 - alpha])
    val_l = med - lower
    val_u = upper - med
    return med, val_l, val_u


def extract_npercent_from_best(samples, percent=68, best=None):
    """

    Returns
    -------
    best : float
    err_lower : float
        best - lower_bound
    err_upper : float
        upper_bound - best
    """
    samples = np.asarray(samples)
    if best is None:
        best = np.median(samples)

    lower_target = percent / 2.0
    upper_target = 100 - percent / 2.0

    cdf = np.sort(samples)
    n = len(cdf)

    frac = np.searchsorted(cdf, best) / n * 100

    if frac >= lower_target and frac <= upper_target:
        lo, hi = np.percentile(samples, [100 - upper_target, upper_target])
    elif frac < lower_target:
        lo = np.min(samples)
        hi = np.percentile(samples, percent)
    else:
        lo = np.percentile(samples, 100 - percent)
        hi = np.max(samples)

    return best, best - lo, hi - best



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
        "--nsigma", type=float, default=1.0,
        help="n-sigma uncertainty")
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

    print(f"  Use equation in {args.paper} with nsigma of {args.nsigma}")

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

    # Add n-sigma
    chi2_nsigma = calc_confidence_chi2(args.paper, chi2_min, dof, args.nsigma, args.reduce)

    print(f"  Extract solutions with chi2 < {chi2_min:.4f} + {chi2_nsigma:.4f} = {chi2_min+chi2_nsigma:.4f}")

    # Uncertainties of TIs 
    chi2_arr = np.array(df["chi2"])
    #   TI of regolith with chi2 < chi2_min + chi2_nsigma
    TIrego_arr = np.array(df["TIrego"])
    TIrego_arr_sig = TIrego_arr[chi2_arr < chi2_min + chi2_nsigma]
    TIrego_nsigl, TIrego_nsigu = np.min(TIrego_arr_sig), np.max(TIrego_arr_sig)
    #   TI of rock with chi2 < chi2_min + chi2_3sigma
    TIrock_arr = np.array(df["TIrock"])
    TIrock_arr_sig = TIrock_arr[chi2_arr < chi2_min + chi2_nsigma]
    TIrock_nsigl, TIrock_nsigu = np.min(TIrock_arr_sig), np.max(TIrock_arr_sig)
    #   Htheta with chi2 < chi2_min + chi2_3sigma
    Htheta_arr = np.array(df["Htheta"])
    Htheta_arr_sig = Htheta_arr[chi2_arr < chi2_min + chi2_nsigma]
    Htheta_nsigl, Htheta_nsigu = np.min(Htheta_arr_sig), np.max(Htheta_arr_sig)
    #   alpha with chi2 < chi2_min + chi2_3sigma
    alpha_arr = np.array(df["alpha"])
    alpha_arr_sig = alpha_arr[chi2_arr < chi2_min + chi2_nsigma]
    alpha_nsigl, alpha_nsigu = np.min(alpha_arr_sig), np.max(alpha_arr_sig)

    text = (
        f"TIrego = ${TIrego_min}_" + "{" + f"-{TIrego_min-TIrego_nsigl}" + "}^" 
        "{" + f"+{TIrego_nsigu-TIrego_min}" + "}" + f"$ (N={len(TIrego_arr_sig)})\n"
        f"TIrock = ${TIrock_min}_" + "{" + f"-{TIrock_min-TIrock_nsigl}" + "}^" 
        "{" + f"+{TIrock_nsigu-TIrock_min}" + "}" + f"$ (N={len(TIrock_arr_sig)})\n"
        f"Htheta = ${Htheta_min}_" + "{" + f"-{Htheta_min-Htheta_nsigl}" + "}^" 
        "{" + f"+{Htheta_nsigu-Htheta_min}" + "}" + f"$ (N={len(Htheta_arr_sig)})\n"
        f"alpha = ${alpha_min}_" + "{" + f"-{alpha_min-alpha_nsigl:.2f}" + "}^" 
        "{" + f"+{alpha_nsigu-alpha_min:.2f}" + "}" + f"$ (N={len(alpha_arr_sig)})\n"
        )

    param_cols = ['Htheta', 'TIrego', 'TIrock', 'alpha']

    # Choose reliable once
    df1 = df[df["chi2"] < chi2_min + chi2_nsigma]

    # Remove specific columns to avoid a following error
    # > ValueError: It looks like the parameter(s) in column(s) 1 have no dynamic range. Please provide a `range` argument.
    for p in param_cols:
        N_p = len(set(df1[p]))
        if N_p == 1:
            param_cols.remove(p)
            print(f"Not plot {p}")


    data_array = df1[param_cols].values
    fig = corner.corner(
        data_array, labels=param_cols, 
        #show_titles=True, 
        label_kwargs={"fontsize": 10}, 
        #title_kwargs={"fontsize": 12}, 
        smoonth=1,
        plot_datapoints=True,     
        #plot_density=True,       
        #plot_contours=True,
        plot_density=False,       
        plot_contours=False,
        bins=50
        )

    percent_interest = 68
    #fig.text(0.55, 0.8, f"The best fit (not median) and the interval\nthat contains {percent_interest}% of the samples are shown.")
    fig.text(0.55, 0.8, f"Note: All data points are\nincluded within the range of figures.")
    
    # Obtain axes 
    axes = np.array(fig.axes).reshape(len(param_cols), len(param_cols))
    # Show median and +- percent/2
    for i, col in enumerate(param_cols):
        ax = axes[i, i]

        if col == "TIrego":
            val_chi2_min = TIrego_min
        elif col == "TIrock":
            val_chi2_min = TIrock_min
        elif col == "Htheta":
            val_chi2_min = Htheta_min
        elif col == "alpha":
            val_chi2_min = alpha_min

        # 1. Calculate median and 1-sigma uncertainties to include percent_interest% samples
        #med, val_l, val_u = extract_npercent(data_array[:, i], percent_interest)
        #print(f"med, val_l, val_u = {med:.1f}, {val_l:.1f}, {val_u:.1f}")

        # 2. de Kleer+2024 use not median but the best fit value
        # Best fit + percent_interest 
        _, val_l, val_u = extract_npercent_from_best(data_array[:, i], percent_interest, best=val_chi2_min)
        #print(f"bestfit, val_l, val_u = {val_chi2_min:.1f}, {val_l:.1f}, {val_u:.1f}")
        #text = f"{col} = ${val_chi2_min:.1f}_" + "{" + f"-{val_l:.1f}" + "}^" + "{" + f"+{val_u:.1f}" + "}$"

        # 3. Use all samples because already extracted with chi-squared
        val_arr = np.array(df[col])
        val_arr_sig = val_arr[chi2_arr < chi2_min + chi2_nsigma]
        val_nsigl, val_nsigu = np.min(val_arr_sig), np.max(val_arr_sig)
        text = (
                f"{col} = ${val_chi2_min:.2f}_" + "{" + f"-{val_chi2_min-val_nsigl:.2f}" + "}^" 
                "{" + f"+{val_nsigu-val_chi2_min:.2f}" + "}$")

        ax.axvline(val_chi2_min, color="red", linestyle="solid", linewidth=1.5, label="Best fit", zorder=100)
        #ax.axvline(val_chi2_min-val_l, color="red", linestyle="dashed", linewidth=2,)
        #ax.axvline(val_chi2_min+val_u, color="red", linestyle="dashed", linewidth=2,)
        ax.set_title(text, fontsize=12)
        ax.legend(fontsize=10)
        # Add margin
        #ax.set_xlim(np.array(ax.get_xlim()) * [0.8, 1.2])


    fig.savefig(args.out, dpi=300, bbox_inches='tight')
