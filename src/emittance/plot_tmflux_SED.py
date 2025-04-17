#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot results of tmflux.
"""
from argparse import ArgumentParser as ap
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt

from myplot import mycolor, mymark


if __name__ == "__main__":
    parser = ap(description="Plot results of tmflux.")
    parser.add_argument(
        "res", type=str, 
        help="Results of do_tmflux_MC.py")
    parser.add_argument(
        "--key_flux", type=str, nargs="*",
        help="Key to specify a flux of interest")
    parser.add_argument(
        "--flux_th", type=float, nargs="*",
        help="Threshold of flux in Jy to divide samples")
    parser.add_argument(
        "--out", type=str, default="tmflux.jpg",
        help="Output file")
    args = parser.parse_args()
   

    df = pd.read_csv(args.res, sep=" ")
    
    # Diameter in m
    df["D"] = df["D"]*1e3
    df0 = df.copy()
    df1 = df.copy()
    
    for (key, th) in zip(args.key_flux, args.flux_th):
        df0 = df0[df0[key] < th]
        df1 = df1[df1[key] >= th]

    # Pick up example
    N_sample = 1
    np.random.seed(0)
    df0 = df0.sample(n=N_sample)
    df1 = df1.sample(n=N_sample)
    df0 = df0.reset_index()
    df1 = df1.reset_index()


    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_axes([0.15, 0.15, 0.7, 0.7])
    ax.set_xlabel("Wavelength [micron]")
    # Flux in mJy
    ax.set_ylabel("Flux density [mJy]")

    col0, mark0 = mycolor[0], mymark[0]
    col1, mark1 = mycolor[1], mymark[1]
    col_th, lw_th = "black", 3
    
    for idx_w, (key, th) in enumerate(zip(args.key_flux, args.flux_th)):

        D0   = df0["D"][0]
        D1   = df1["D"][0]
        pV0   = df0["pv"][0]
        pV1   = df1["pv"][0]
        eta0 = df0["eta"][0]
        eta1 = df1["eta"][0]

        # Plot label once
        if idx_w == 0:
            label0 = (
                f"Faint example\n  "
                f"D={D0:.1f} m\n  "
                f"$p_V$={pV0:.2f} \n  "
                f"$\eta$={eta0:.1f}"
                )
            label1 = (
                f"Bright example\n  "
                f"D={D1:.1f} m\n  "
                f"$p_V$={pV1:.2f} \n  "
                f"$\eta$={eta1:.1f}"
                )
        else:
            label0 = None
            label1 = None

        # Extract wavelength
        # flux8 -> 8,7
        if key == "flux8":
            w = 8.7
        elif key == "flux10":
            w = 10.7
        elif key == "flux11":
            w = 11.7
        df0["w"] = w
        df1["w"] = w
        
        # Model examples
        ax.scatter(
            df0["w"], df0[key]*1e3, color=col0, marker=mark0, facecolor="None", label=label0)
        ax.scatter(df1["w"], df1[key]*1e3, color=col1, marker=mark1, facecolor="None", label=label1)

        # Upper limit
        start, end = (w, th*1000), (w, th*1000-3)
        ax.annotate(
            '', xy=end, xytext=start,size=10,
            arrowprops=dict(arrowstyle='-|>',
            connectionstyle='arc3', lw=lw_th, facecolor=col_th, edgecolor=col_th)
        )
        w0, w1 = w-0.1, w+0.1 
        ax.hlines(th*1000, w0, w1, lw=lw_th, color=col_th)
        
    ax.set_xlim([8, 12])
    ax.legend(ncol=2).get_frame().set_alpha(1.0)
    plt.savefig(args.out)
