#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plot histograms of estimated parameters.
"""
from argparse import ArgumentParser as ap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from trem.emittance.common_emittance import calc_confidence_chi2


def plot_hist(df_list, key, label_list, l=0.25, u=0.75, bins=30, out="hist.jpg"):
    fig = plt.figure()
    ax = fig.add_axes([0.20, 0.15, 0.75, 0.8])
    
    for idx, df in enumerate(df_list):
        if len(df) > 0:
            lv, uv = np.quantile(df[key], [l, u])
            label = f"{label_list[idx]} {lv:.2f}--{uv:.2f} ({l*100}--{u*100}\%) N={len(df)}"
            ax.hist(df[key], histtype="step", density=True, bins=bins, label=label)
    ax.legend()
    ax.set_xlabel(key)
    ax.set_yscale("log")
    ax.set_ylabel("Normalized counts")

    plt.savefig(out)


if __name__ == "__main__":
    parser = ap(description="Plot histogram of parameters")
    parser.add_argument(
        "res", type=str,
        help="Result of TPM with chi2.")
    parser.add_argument(
        "-x", type=str, default="TIrock", 
        help="x axis of the plot (A, TIrock, TIrego Htheta, etc.)")
    parser.add_argument(
        "--dof", type=int, default=1,
        help="Degree of freedom")
    parser.add_argument(
        "--out", type=str, default="hist.jpg",
        help="Output file")
    args = parser.parse_args()
   
    
    # Read files 
    df = pd.read_csv(args.res, sep=" ")
    dof = args.dof
    
    # reduce
    df["chi2"] = df["chi2"]/dof
    chi2_min = np.min(df["chi2"])
    print(f"reduced chi2 minimum = {chi2_min:.3f}")


    ## 1. Used in Polishool 2014 etc.
    chi2_1sigma_P14 = calc_confidence_chi2("P14", chi2_min, dof, 1, True)
    chi2_3sigma_P14 = calc_confidence_chi2("P14", chi2_min, dof, 3, True)
    print(f"P14 (1sigma, 3sigma) = ({chi2_1sigma_P14:.2f}, {chi2_3sigma_P14:.2f})")
    
    ## 2. Used in Vokrouhlicky+2017 etc.
    chi2_1sigma_V17 = calc_confidence_chi2("V17", chi2_min, dof, 1, True)
    chi2_3sigma_V17 = calc_confidence_chi2("V17", chi2_min, dof, 3, True)
    print(f"V17 (1sigma, 3sigma) = ({chi2_1sigma_V17:.2f}, {chi2_3sigma_V17:.2f})")


    # Extract samples for Case 1.
    ## 1 sigma
    df_P14_1 = df[df["chi2"] < chi2_min + chi2_1sigma_P14]
    print(f"P14, 1sigma, N={len(df_P14_1)}")
    ## 3 sigma
    df_P14_3 = df[df["chi2"] < chi2_min + chi2_3sigma_P14]
    print(f"P14, 3sigma, N={len(df_P14_3)}")
    
    # Extract samples for Case 2.
    ## 1 sigma
    df_V17_1 = df[df["chi2"] < chi2_min + chi2_1sigma_V17]
    print(f"V17, 1sigma, N={len(df_V17_1)}")
    ## 3 sigma
    df_V17_3 = df[df["chi2"] < chi2_min + chi2_3sigma_V17]
    print(f"V17, 3sigma, N={len(df_V17_3)}")
    
    # All solution
    plot_hist([df_P14_1, df_P14_3, df_V17_1, df_V17_3], args.x, ["P14, 1sigma", "P14, 3sigma", "V17, 1sigma", "V17, 3 sigma"], l=0.25, u=0.75, bins=30, out=args.out)
