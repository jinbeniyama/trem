#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plot results of tmflux.
"""
from argparse import ArgumentParser as ap
import numpy as np
import pandas as pd 
import scipy.stats as stats
import matplotlib.pyplot as plt
from matplotlib.patches import ArrowStyle

from myplot import mycolor, myls

def get_range(arr, percentage):
    # Median, 1-sigma
    v_med = np.median(arr)
    v_l = stats.scoreatpercentile(arr, (100.-percentage)/2.)
    v_u = stats.scoreatpercentile(arr, percentage+(100.-percentage)/2.)
    return v_med, v_l, v_u


def plot_mcres(f, key_flux, flux_th, nbin=50, out=None):
    df = pd.read_csv(f, sep=" ")

    df["D"] = df["D"]*1e3
    df0 = df[df[key_flux] < flux_th]
    df1 = df[df[key_flux] >= flux_th]

    Nall, N0, N1 = len(df), len(df0), len(df1)

    fig = plt.figure(figsize=(16, 12))
    # Upper: inputs
    # Lower: outputs
    xwi, ywi = 0.22, 0.22
    colall = "black"
    col0, col1 = mycolor[0], mycolor[1]
    ls0, ls1 = myls[0], myls[1]
    mark0, mark1 = ".", "+"
    s0, s1 = 5, 5

    lw_range = 2
    ls_range = "dashed"


    # input H
    ax_u1 = fig.add_axes([0.10, 0.70, xwi, ywi])
    ax_u1.set_xlabel("input H")
    ax_u1.set_ylabel("N")
    ax_u1.hist(df["H"], histtype="step", bins=nbin, rwidth=0.6, color=colall, label=f"N={Nall}")
    # input pv
    ax_u2 = fig.add_axes([0.40, 0.70, xwi, ywi])
    ax_u2.set_xlabel("input pv")
    ax_u2.set_ylabel("N")
    ax_u2.hist(df["pv"], histtype="step", bins=nbin, rwidth=0.6, color=colall, label=f"N={Nall}")
    # input eta
    ax_u3 = fig.add_axes([0.70, 0.70, xwi, ywi])
    ax_u3.set_xlabel("input eta")
    ax_u3.set_ylabel("N")
    ax_u3.hist(df["eta"], histtype="step", bins=nbin, rwidth=0.6, color=colall, label=f"N={Nall}")



    percentage = 68.2
    # D
    ax_m1 = fig.add_axes([0.10, 0.40, xwi, ywi])
    ax_m1.set_xlabel("D [m]")
    ax_m1.set_ylabel("N")
    ax_m1.hist(df0["D"], histtype="step", bins=nbin, color=col0, ls=ls0)
    ax_m1.hist(df1["D"], histtype="step", bins=nbin, color=col1, ls=ls1)
    ymax_d  = ax_m1.get_ylim()[1]
    # Median, 1-sigma
    d_median, d_l, d_u = get_range(df0["D"], percentage)
    derr_l = d_median - d_l
    derr_u = d_u - d_median
    print(f"Diameter: (median, +1sigma, -1sigma) = ({d_median:.1f}, {derr_u:.1f}, {derr_l:.1f})\n")
    ax_m1.hlines(ymax_d/3, d_l, d_u, color=col0, lw=lw_range, ls=ls0)
    start, end = (d_l, ymax_d/3), (d_u, ymax_d/3)
    ax_m1.annotate(
        '', xy=end, xytext=start, size=20,
        arrowprops=dict(arrowstyle=ArrowStyle('|-|', widthA=0.5, widthB=0.5),
        shrinkA=0, shrinkB=0, lw=lw_range, ls=ls0,
        connectionstyle='arc3', facecolor=col0, color=col0)
        )
    label0  = f"flux < {flux_th*1000} mJy" + "\n" + f"D={d_median:.1f}" + "$_{" + f"-{derr_l:.1f}" + "}^{" + f"+{derr_u:.1f}" + "} (1\sigma)$"
    ax_m1.scatter(
        d_median, ymax_d/3, color=col0, marker="o", s=50, zorder=100, label=label0)

    # pv
    ax_m2 = fig.add_axes([0.40, 0.40, xwi, ywi])
    ax_m2.set_xlabel("pv")
    ax_m2.set_ylabel("N")
    ax_m2.hist(df0["pv"], histtype="step", bins=nbin, color=col0, ls=ls0)
    ax_m2.hist(df1["pv"], histtype="step", bins=nbin, color=col1, ls=ls1)
    ymax_pv  = ax_m2.get_ylim()[1]
    pv_median, pv_l, pv_u = get_range(df0["pv"], percentage)
    pverr_l = pv_median - pv_l
    pverr_u = pv_u - pv_median
    print(f"pv: (median, +1sigma, -1sigma) = ({pv_median:.1f}, {pverr_u:.1f}, {pverr_l:.1f})\n")
    ax_m2.hlines(ymax_pv/3, pv_l, pv_u, color=col0, lw=lw_range, ls=ls0)
    start, end = (pv_l, ymax_pv/3), (pv_u, ymax_pv/3)
    ax_m2.annotate(
        '', xy=end, xytext=start, size=20,
        arrowprops=dict(arrowstyle=ArrowStyle('|-|', widthA=0.5, widthB=0.5),
        shrinkA=0, shrinkB=0, lw=lw_range, ls=ls0,
        connectionstyle='arc3', facecolor=col0, color=col0)
        )
    label0 = f"flux < {flux_th*1000} mJy" + "\n" + f"p$_V$={pv_median:.2f}" + "$_{" + f"-{pverr_l:.2f}" + "}^{" + f"+{pverr_u:.2f}" + "} (1\sigma)$"
    ax_m2.scatter(
        pv_median, ymax_pv/3, color=col0, marker="o", s=50, zorder=100, label=label0)

    # eta
    ax_m3 = fig.add_axes([0.70, 0.40, xwi, ywi])
    ax_m3.set_xlabel("eta")
    ax_m3.set_ylabel("N")
    ax_m3.hist(df0["eta"], histtype="step", bins=nbin, color=col0, ls=ls0)
    ax_m3.hist(df1["eta"], histtype="step", bins=nbin, color=col1, ls=ls1)
    ymax_eta  = ax_m3.get_ylim()[1]
    eta_median, eta_l, eta_u = get_range(df0["eta"], percentage)
    etaerr_l = eta_median - eta_l
    etaerr_u = eta_u - eta_median
    print(f"eta: (median, +1sigma, -1sigma) = ({eta_median:.1f}, {etaerr_u:.1f}, {etaerr_l:.1f})\n")
    ax_m3.hlines(ymax_eta/3, eta_l, eta_u, color=col0, lw=lw_range, ls=ls0)
    start, end = (eta_l, ymax_eta/3), (eta_u, ymax_eta/3)
    ax_m3.annotate(
        '', xy=end, xytext=start, size=20,
        arrowprops=dict(arrowstyle=ArrowStyle('|-|', widthA=0.5, widthB=0.5),
        shrinkA=0, shrinkB=0, lw=lw_range, ls=ls0,
        connectionstyle='arc3', facecolor=col0, color=col0)
        )
    label0 = f"flux < {flux_th*1000} mJy" + "\n" +f"$\eta$={eta_median:.2f}" + "$_{" + f"-{etaerr_l:.2f}" + "}^{" + f"+{etaerr_u:.2f}" + "} (1\sigma)$"
    ax_m3.scatter(
        eta_median, ymax_eta/3, color=col0, marker="o", s=50, zorder=100, label=label0)


    # D vs. pv
    ax_l1 = fig.add_axes([0.10, 0.1, xwi, ywi])
    ax_l1.set_xlabel("D [m]")
    ax_l1.set_ylabel("pv")
    ax_l1.scatter(df0["D"], df0["pv"], color=col0, marker=mark0, s=s0, label=f"N={N0}")
    ax_l1.scatter(df1["D"], df1["pv"], color=col1, marker=mark1, s=s1, label=f"N={N1}")
    ax_l1.set_xscale("log")
    ax_l1.set_yscale("log")

    # eta vs. pv
    ax_l2 = fig.add_axes([0.40, 0.1, xwi, ywi])
    ax_l2.set_xlabel("eta")
    ax_l2.set_ylabel("pv")
    ax_l2.scatter(df0["eta"], df0["pv"], color=col0, marker=mark0, s=s0, label=f"N={N0}")
    ax_l2.scatter(df1["eta"], df1["pv"], color=col1, marker=mark1, s=s1, label=f"N={N1}")

    # D vs. eta
    ax_l3 = fig.add_axes([0.70, 0.1, xwi, ywi])
    ax_l3.set_xlabel("D [m]")
    ax_l3.set_ylabel("eta")
    ax_l3.scatter(df0["D"], df0["eta"], color=col0, marker=mark0, s=s0, label=f"N={N0}")
    ax_l3.scatter(df1["D"], df1["eta"], color=col1, marker=mark1, s=s1, label=f"N={N1}")

    for ax in fig.axes:
        ax.legend()
    if out:
        plt.savefig(out)


if __name__ == "__main__":
    parser = ap(description="Plot results of tmflux.")
    parser.add_argument(
        "res", type=str, 
        help="Results of do_tmflux_MC.py")
    parser.add_argument(
        "--key_flux", type=str, default="flux8",
        help="Key to specify a flux of interest")
    parser.add_argument(
        "--flux_th", type=float, default=None,
        help="Threshold of flux in Jy to divide samples")
    parser.add_argument(
        "--nbin", type=int, default=50,
        help="Number of bins of histograms")
    parser.add_argument(
        "--out", type=str, default="tmflux.jpg",
        help="Output file")
    args = parser.parse_args()
   
    plot_mcres(
        args.res, args.key_flux, flux_th=args.flux_th, nbin=args.nbin, out=args.out)
