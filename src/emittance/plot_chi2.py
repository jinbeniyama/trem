#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot results of TPM.

w/ fixscale argument
    Can fix scale parameter (recalculate chi-squared in the code)

w/ scale_per_obs argument
    Can introduce new scale factors 's_j' (j: observation time) in this code. 
    No scale factor is introduced for a single observations (i.e., photometry).

Example
-------
# A vs. chi square
> plot_tpm.py tpmout* -x A --out result_A.jpg
# TI (thermal inertia) vs. chi square
> plot_tpm.py tpmout* -x TI --out result_TI.jpg
"""
import os 
from argparse import ArgumentParser as ap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from myplot import mycolor, mymark

from trem.common import mycolor, mymark
from trem.emittance.common_emittance import extract_flux, calc_chi2, introduce_var_scalefactor, calc_confidence_chi2


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
        "--paper", type=str, default="P14",
        help="P14 or V17, type of uncertainty")
    parser.add_argument(
        "--logx", action="store_true", default=False,
        help="Use log scale in x axis")
    parser.add_argument(
        "--logy", action="store_true", default=False,
        help="Use log scale in y axis")
    parser.add_argument(
        "--out", type=str, default="result.jpg",
        help="Output file")
    parser.add_argument(
        "--out_df", type=str, default=None,
        help="Output chi2 etc as dataframe")
    parser.add_argument(
        "--outdir", type=str, default=".",
        help="Directory for output file")
    args = parser.parse_args()
   
    
    # Parse arguments =========================================================
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    fixscale = args.fixscale
    # Parse arguments =========================================================


    # Read files 
    df = pd.read_csv(args.res, sep=" ")
    # Add dummy
    df["scalefactor"] = 1
    Htheta_list_sort = sorted(list(set(df["Htheta"])))
    TI_list_sort = sorted(list(set(df["TI"])))
    df_temp = df[
        (df["Htheta"] == Htheta_list_sort[0]) & 
        (df["TI"] == TI_list_sort[0])
        ]
    # Number of data points
    N_data = len(df_temp)
    
    # Parameter of interest in X label
    val = args.x
    if val == "A":
        xlabel = "Bond albedo"
    elif val == "TI":
        xlabel = "Thermal intertia [tiu]"
    elif val == "CA":
        xlabel = "Crater opening angle [deg]"
    elif val == "CR":
        xlabel = "Crater covering ratio"
     

    # Calculate dof ===========================================================
    # chi2_reduce = chi2/dof, where dof is a degree of freedom
    # dof = N_data - N_param, where N_data and N_param are numbers of 
    # data points and parameters, respectively
    # In Marco's code, fit parameter is only diameter (i.e., N_param = 1)
    # Thus, dof = N_data - 1

    # Calculate degree of freedom
    N_param = args.N_param
    dof = N_data - N_param
    print(f"  Calculate reduced chi square with dof = {dof}")
    # Calculate dof ===========================================================
    

    # Plot
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_axes([0.15, 0.15, 0.55, 0.80])
    ax.grid(which="both", color="gray",linewidth=0.2)

    axin = ax.inset_axes([0.60, 0.50, 0.35, 0.35])

    if args.logx:
        ax.set_xscale("log")
    if args.logy:
        ax.set_yscale("log")

    # Use reduced chi2
    if args.reduce:
        ylabel = r"Reduced $\chi_2$  ($\nu$=" + f"{N_data}-{N_param}={dof})"
    else:
        ylabel = r"$\chi_2$"

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if args.ylim:
        y0, y1 = args.ylim
        ax.set_ylim([y0, y1])

    # Calculate chi2 fixing scale factor to 1 =================================
    if fixscale:
        chi2_new_list = []
        TI_new_list, Htheta_new_list = [], []
        Htheta_used = []
        for idx_Htheta, Htheta in enumerate(Htheta_list_sort):
            print(f"idx_Htheta = {idx_Htheta+1}/{len(Htheta_list_sort)}")
            for idx_TI, TI in enumerate(TI_list_sort):

                df_temp = df[
                    (df["Htheta"] == Htheta) &
                    (df["TI"] == TI)] 
                chi2 = calc_chi2(df_temp)
                # Use reduced chi2
                if args.reduce:
                    chi2 = chi2/dof 
                else:
                    pass

                # Index of Htheta 
                idx_Htheta = Htheta_list_sort.index(Htheta)
                col, mark = mycolor[idx_Htheta], mymark[idx_Htheta]
                if Htheta not in Htheta_used:
                    label = f"Hapke angle = ({Htheta:.1f})"
                    Htheta_used.append(Htheta)
                else:
                    label = None
                for a in [ax, axin]:
                    a.scatter(
                        np.min(df_temp[val]), chi2, color=col, marker=mark, s=50, 
                        facecolor="None", label=label)
                chi2_new_list.append(chi2)
                TI_new_list.append(TI)
                Htheta_new_list.append(Htheta)

        # Add global minimum chi2
        chi2_min = np.min(chi2_new_list)
        # Index of global chi2_min
        idx_min = chi2_new_list.index(min(chi2_new_list))
        chi2_arr = np.array(chi2_new_list)
    # Calculate chi2 fixing scale factor to 1 =================================


    # Calculate chi2 with scale factors per observation =======================
    elif args.scale_per_obs:
        assert False, "Not yet implemented"
        chi2_new_list = []
        TI_new_list, Htheta_new_list = [], []
        Htheta_list_sort = sorted(list(set(Htheta_list)))
        Htheta_used = []
        for f in args.res:
            # To check the origin of genmesh FAILED
            #print(f)
            # Extract fluxes with scale factor = 1
            df_temp = extract_flux(f, fixscale=True)
            # Introduce variable scale factors
            df_temp = introduce_var_scalefactor(df_temp)
            
            # Then calculate chi2 with 'f_obs', 'scalefactor', 'f_model', 'ferr_obs'.
            # 'scele factor' is not a constant any more
            sf_list = set(df_temp["scalefactor"])
            N_sf = len(sf_list)
            #print(f"  Number of scale parameters: N_sf = {N_sf}")
            #print(f"     {sf_list}")
            chi2 = calc_chi2(df_temp)
            # Use reduced chi2
            if args.reduce:
                chi2 = chi2/dof 
            else:
                pass

            Htheta = np.min(df_temp["Htheta"])
            # Index of Htheta 
            idx_Htheta = Htheta_list_sort.index(Htheta)
            col, mark = mycolor[idx_Htheta], mymark[idx_Htheta]
            if Htheta not in Htheta_used:
                label = f"Hapke angle = ({Htheta:.1f})"
                Htheta_used.append(Htheta)
            else:
                label = None
            for a in [ax, axin]:
                a.scatter(
                    np.min(df_temp[val]), chi2, color=col, marker=mark, s=50, 
                    facecolor="None", label=label)
            chi2_new_list.append(chi2)
            TI_new_list.append(TI)
            Htheta_new_list.append(Htheta)

        # Add global minimum chi2
        chi2_min = np.min(chi2_new_list)
        # Index of global chi2_min
        idx_min = chi2_new_list.index(min(chi2_new_list))
        chi2_arr = np.array(chi2_new_list)

    # Save them
    if args.out_df:
        df_out = pd.DataFrame({
            "chi2": chi2_new_list,
            "TI": TI_new_list,
            "Htheta": Htheta_new_list,
            })
        out_df = os.path.join(args.outdir, args.out_df)
        df_out.to_csv(out_df, sep=" ")

    # Calculate chi2 with scale factors per observation =======================
    

    # Add 1-sigma, 3-sigma
    chi2_1sigma = calc_confidence_chi2(args.paper, chi2_min, dof, 1, args.reduce)
    chi2_3sigma = calc_confidence_chi2(args.paper, chi2_min, dof, 3, args.reduce)
    
    for a in [ax, axin]:
        xmin, xmax = a.get_xlim()
        a.hlines(chi2_min + chi2_1sigma, xmin, xmax, ls="dashed", color="black", label=r"1$\sigma$" + f" ({chi2_1sigma:.2f}) {args.paper}")
        a.hlines(chi2_min + chi2_3sigma, xmin, xmax, ls="dotted", color="black", label=r"3$\sigma$" + f" ({chi2_3sigma:.2f}) {args.paper}")
        a.set_xlim([xmin, xmax])

    ## Extract parameters which give chi2_min
    TI_min = TI_new_list[idx_min]
    Htheta_min = Htheta_new_list[idx_min]
    text = f"Minimum chi2 = {chi2_min:.2f} w/ (TI, Htheta) = ({TI_min:.2f}, {Htheta_min:.2f})"
    ax.text(0.1, 0.95, text, size=12, transform=ax.transAxes)
    ax.legend(bbox_to_anchor=(1.05, 1), borderaxespad=0, fontsize=8, ncol=2)

    # Small axis
    # Uncertainties of value of interest
    if val == "TI":
        val_arr = np.array(TI_new_list)
        valbest = TI_min
    else:
        assert False, "Not implimented"

    val_arr_sig = val_arr[chi2_arr < chi2_min + chi2_3sigma]
    val3sigl, val3sigu = np.min(val_arr_sig), np.max(val_arr_sig)
    text = (
        f"{val}= ${valbest:.2f}_" + "{" + f"-{valbest-val3sigl:.2f}" + "}^" 
        "{" + f"+{val3sigu-valbest:.2f}" + "}" + f"$ (N={len(val_arr_sig)})"
        )
    axin.text(0.1, 0.80, text, size=12, transform=axin.transAxes)
    val_0, val_1 = val3sigl*0.8, val3sigu*1.2
    chi2_0, chi2_1 = chi2_min*0.8, (chi2_min + chi2_3sigma)*1.2
    axin.set_xlim([val_0, val_1])
    axin.set_ylim([chi2_0, chi2_1])

    out = os.path.join(args.outdir, args.out)
    fig.savefig(out)
