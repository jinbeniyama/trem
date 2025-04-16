#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plot results of TPM using a brute-force method.

w/ fixscale argument
    Can fix scale parameter (recalculate chi-squared in the code)

w/ scale_per_obs argument
    Can introduce new scale factors 's_j' (j: observation time) in this code. 
    No scale factor is introduced for a single observations (i.e., photometry).

Example
-------
# A vs. chi square
> plot_tpm_brute.py tpmout* -x A --out result_A.jpg
# TI (thermal inertia) vs. chi square
> plot_tpm_brute.py tpmout* -x TI --out result_TI.jpg
"""
import os 
from argparse import ArgumentParser as ap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from myplot import mycolor, mymark

from tpmwrapper.common import extract_flux, calc_chi2, introduce_var_scalefactor, calc_confidence_chi2


if __name__ == "__main__":
    parser = ap(description="Plot results of TPM.")
    parser.add_argument(
        "res", nargs="*", 
        help="Results of rumtpm")
    parser.add_argument(
        "-x", type=str, default="TI", 
        help="x axis of the plot (A, TI)")
    parser.add_argument(
        "--ylim", nargs=2, type=float, default=None, 
        help="x axis of the plot")
    parser.add_argument(
        "--fixscale", action="store_true", default=False,
        help="Fix scale factor to 1.")
    parser.add_argument(
        "--scale_per_obs", action="store_true", default=False,
        help="Use scale factors per observation")
    parser.add_argument(
        "--reduce", action="store_true", default=False,
        help="Reduced chi square")
    parser.add_argument(
        "--N_param", type=int, default=1,
        help="Number of parameters")
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
        "--outdir", type=str, default=".",
        help="Directory for output file")
    args = parser.parse_args()
   
    
    # Parse arguments =========================================================
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    fixscale = args.fixscale
    # Parse arguments =========================================================
    
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
     

    # Read files and extract physical parameters and chi2 =====================
    A_list, TI_list, CA_list, CR_list = [], [], [], []
    scale_list, Dv_list, Ds_list = [], [], []
    Htheta_list = []
    chi2_list = []
    N_data = 0
    
    # Note: All TPM results should be from the same conditions 
    #       (obs file, ephem file) with different parameters (A, TI, etc.).
    for idx_res, res in enumerate(args.res):
        with open (res, "r") as f:
            print(f"Read {res}")
            lines = f.readlines()

            # Extract line starts with r> (idx of -6) 
            for l in lines:
                if l[:2] == "r>":
                    l_split = l.split()
                    A = l_split[5]
                    A_list.append(float(A))
                    TI = l_split[1]
                    TI_list.append(float(TI))
                    CA = l_split[3]
                    CA_list.append(float(CA))
                    CR = l_split[4]
                    CR_list.append(float(CR))
                    Htheta = l_split[2]
                    Htheta_list.append(float(Htheta))
                    # scale factor
                    scale = l_split[6]
                    scale_list.append(float(scale))
                    # Dv: ??
                    Dv = l_split[7]
                    Dv_list.append(float(Dv))
                    # Ds: ??
                    Ds = l_split[8]
                    Ds_list.append(float(Ds))
                    # chi2 (not reduced yet)
                    chi2 = l_split[10]
                    chi2_list.append(float(chi2))
                elif (idx_res == 0) and ("ReadObs: Observation" in l):
                    # This selection is valid for both TI==0 and TI!=0
                    # Extract Ndata
                    l_split = l.split()
                    N_obs = int(l_split[4].split("=")[1])
                    N_data += N_obs

                    
            print(f"    Physical parameters extracted.")

    print("")
    print(f"Number of observations is {N_data}")
    print("")
    # Count nan
    nan_A    = np.isnan(A_list).sum()
    nan_TI   = np.isnan(TI_list).sum()
    nan_CA   = np.isnan(CA_list).sum()
    nan_CR   = np.isnan(CR_list).sum()
    nan_sc   = np.isnan(scale_list).sum()
    nan_Dv   = np.isnan(Dv_list).sum()
    nan_Ds   = np.isnan(Ds_list).sum()
    nan_chi2 = np.isnan(chi2_list).sum()

    print("Number of nan.")
    # These are fixed parameters. (i.e., always wo/ nan)
    #print(f"    nan in A    : {nan_A}")
    #print(f"    nan in TI   : {nan_TI}")
    #print(f"    nan in CA   : {nan_CA}")
    #print(f"    nan in CR   : {nan_CR}")
    # N of nan is always the same?
    print(f"    nan in sc   : {nan_sc}")
    print(f"    nan in Dv   : {nan_Dv}")
    print(f"    nan in Ds   : {nan_Ds}")
    print(f"    nan in chi2 : {nan_chi2}")
    # Read files and extract physical parameters and chi2 =====================


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
        Htheta_list_sort = sorted(list(set(Htheta_list)))
        Htheta_used = []
        for f in args.res:
            # To check the origin of genmesh FAILED
            #print(f)
            df_temp = extract_flux(f, fixscale)
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
                label = f"Roughness params.\n  Hapke angle = ({Htheta})"
                Htheta_used.append(Htheta)
            else:
                label = None
            for a in [ax, axin]:
                a.scatter(
                    np.min(df_temp[val]), chi2, color=col, marker=mark, s=50, 
                    facecolor="None", label=label)
            chi2_new_list.append(chi2)

        # Add global minimum chi2
        chi2_min = np.min(chi2_new_list)
        # Index of global chi2_min
        idx_min = chi2_new_list.index(min(chi2_new_list))
        chi2_arr = np.array(chi2_new_list)
    # Calculate chi2 fixing scale factor to 1 =================================


    # Calculate chi2 with scale factors per observation =======================
    elif args.scale_per_obs:
        chi2_new_list = []
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
                label = f"Roughness params.\n  Hapke angle = ({Htheta})"
                Htheta_used.append(Htheta)
            else:
                label = None
            for a in [ax, axin]:
                a.scatter(
                    np.min(df_temp[val]), chi2, color=col, marker=mark, s=50, 
                    facecolor="None", label=label)
            chi2_new_list.append(chi2)

        # Add global minimum chi2
        chi2_min = np.min(chi2_new_list)
        # Index of global chi2_min
        idx_min = chi2_new_list.index(min(chi2_new_list))
        chi2_arr = np.array(chi2_new_list)
    # Calculate chi2 with scale factors per observation =======================
    

    # Use raw chi2 in output files of TPM =====================================
    else:
        if args.reduce:
            chi2_list = [x/dof for x in chi2_list]
        else:
            pass
        # Make a data frame
        df = pd.DataFrame({
            "A": A_list, 
            "TI": TI_list, 
            "CA": CA_list, 
            "CR": CR_list, 
            "chi2": chi2_list}
            )
        # Parameter about roughness is characterized by
        # CA, not CR (i.e., duplicated CRs are used.)
        CAs = df.CA.unique()

        for idx, CA in enumerate(CAs):
            df_temp = df[df["CA"] == CA]
            assert len(df_temp["CA"].unique()) == 1
            assert len(df_temp["CR"].unique()) == 1
            ca, cr = np.min(df_temp["CA"]), np.min(df_temp["CR"])
            col, mark = mycolor[idx], mymark[idx]
            label = f"Roughness params.\n  (CA, CR) = ({ca}, {cr})"
            for a in [ax, axin]:
                a.scatter(
                    df_temp[val], df_temp["chi2"], color=col, marker=mark, s=50, 
                    facecolor="None", label=label)

        # Add global minimum chi2
        chi2_min = np.min(df["chi2"])
        # Index of global chi2_min
        idx_min = chi2_list.index(min(chi2_list))
        chi2_arr = np.array(chi2_list)
    # Use raw chi2 in output files of TPM =====================================


    # Add 1-sigma, 3-sigma
    chi2_1sigma = calc_confidence_chi2(args.paper, chi2_min, dof, 1, args.reduce)
    chi2_3sigma = calc_confidence_chi2(args.paper, chi2_min, dof, 3, args.reduce)
    
    for a in [ax, axin]:
        xmin, xmax = a.get_xlim()
        a.hlines(chi2_min + chi2_1sigma, xmin, xmax, ls="dashed", color="black", label=r"1$\sigma$" + f" ({chi2_1sigma:.2f}) {args.paper}")
        a.hlines(chi2_min + chi2_3sigma, xmin, xmax, ls="dotted", color="black", label=r"3$\sigma$" + f" ({chi2_3sigma:.2f}) {args.paper}")
        a.set_xlim([xmin, xmax])

    # Extract parameters which give chi2_min
    TI_min = TI_list[idx_min]
    A_min = A_list[idx_min]
    Htheta_min = Htheta_list[idx_min]
    CA_min = CA_list[idx_min]
    CR_min = CR_list[idx_min]
    text = f"Minimum chi2 = {chi2_min:.2f} w/ (TI, A, Htheta, CA, CR) = ({TI_min}, {A_min}, {Htheta_min}, {CA_min}, {CR_min})"
    ax.text(0.1, 0.95, text, size=12, transform=ax.transAxes)
    ax.legend(bbox_to_anchor=(1.05, 1), borderaxespad=0)

    # Small axis
    # Uncertainties of value of interest
    if val == "TI":
        val_arr = np.array(TI_list)
        valbest = TI_min
    else:
        assert False, "Not implimented"

    val_arr_sig = val_arr[chi2_arr < chi2_min + chi2_3sigma]
    val3sigl, val3sigu = np.min(val_arr_sig), np.max(val_arr_sig)
    text = (
        f"{val}= ${valbest}_" + "{" + f"-{valbest-val3sigl}" + "}^" 
        "{" + f"+{val3sigu-valbest}" + "}" + f"$ (N={len(val_arr_sig)})"
        )
    axin.text(0.1, 0.80, text, size=12, transform=axin.transAxes)
    val_0, val_1 = val3sigl*0.8, val3sigu*1.2
    chi2_0, chi2_1 = chi2_min*0.8, (chi2_min + chi2_3sigma)*1.2
    axin.set_xlim([val_0, val_1])
    axin.set_ylim([chi2_0, chi2_1])


    out = os.path.join(args.outdir, args.out)
    fig.savefig(out)
