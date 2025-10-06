#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plot results of TPM.

A 2D plot is easy to look at, 
but difficult to interpret because some points are hidden behind others.
"""
from argparse import ArgumentParser as ap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plot_chi2(df, key_x, key_y, key_z, dof, vmin=1, vmax=100, logx=False, logy=True, out=None):
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_axes([0.15, 0.15, 0.70, 0.8])

    ax.set_xlabel(key_x)
    ax.set_ylabel(key_y)
    if logx:
        ax.set_xscale("log")
    if logy:
        ax.set_yscale("log")

    x_list = sorted(list(set(df[key_x])))
    y_list = sorted(list(set(df[key_y])))

    chi2_list = df[key_z]/dof
    sc = ax.scatter(df[key_x], df[key_y], c=chi2_list, vmin=vmin, vmax=vmax)

    chi2_min = np.min(chi2_list)
    ax.text(0.05, 0.05, f"Chi2 minimum: {chi2_min:.2f}", size=20, transform=ax.transAxes)

    plt.colorbar(sc)
    if out:
        plt.savefig(out)


if __name__ == "__main__":
    parser = ap(description="Plot results of TPM.")
    parser.add_argument(
        "res", type=str,
        help="Chi2 etc.")
    parser.add_argument(
        "-x", type=str, default="Htheta", 
        help="x axis of the plot (A, TI)")
    parser.add_argument(
        "-y", type=str, default="TI", 
        help="y axis of the plot (A, TI)")
    parser.add_argument(
        "--reduce", action="store_true", default=False,
        help="Reduced chi square")
    parser.add_argument(
        "--dof", type=int, default=1,
        help="Degree of freedom")
    parser.add_argument(
        "--vr", type=float, nargs=2, default=None,
        help="Value range")
    parser.add_argument(
        "--out", type=str, default="chi2_xy.jpg",
        help="Output file")
    args = parser.parse_args()
   
    
    # Read files 
    df = pd.read_csv(args.res, sep=" ")
    vmin, vmax = args.vr
    plot_chi2(df, args.x, args.y, "chi2", args.dof, vmin=vmin, vmax=vmax, out=args.out)
