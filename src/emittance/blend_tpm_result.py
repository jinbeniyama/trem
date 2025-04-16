#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Search best regolith abundance etc. using dual-component & brute-force method
with fixed cutoff thermal inertia (TI_th).
Introduce a new parameter 'alpha' in this code.

Example
-------
# Calculate best chi2 values with variety of parameters.
> blend_tpm_result.py tpmout* 
"""
import os 
from argparse import ArgumentParser as ap
import numpy as np
import pandas as pd

from tpmwrapper.common import extract_flux, search_regolith_abundance, crater2Htheta


if __name__ == "__main__":
    parser = ap(description="Blend thermal fluxes introducing alpha, regolith abundance")
    parser.add_argument(
        "res", nargs="*", 
        help="Results of rumtpm")
    parser.add_argument(
        "--resdir", type=str, default=".",
        help="Directory in which results are saved")
    parser.add_argument(
        "--TI_thresh", type=float, default=1500,
        help="Threshold of TI of regolith and TI of rocks")
    parser.add_argument(
        "--chi2_min0", type=float, default=200000,
        help="Initial minimum chi2 used to find the best alpha")
    parser.add_argument(
        "--amin", type=float, default=0.0,
        help="Minimum regolith abundance")
    parser.add_argument(
        "--amax", type=float, default=1.0,
        help="Maximum regolith abundance")
    parser.add_argument(
        "--astep", type=float, default=0.1,
        help="Step of regolith abundance")
    parser.add_argument(
        "--fixscale", action="store_true", default=False,
        help="Fix scale factor to 1.")
    parser.add_argument(
        "--fitalpha", action="store_true", default=False,
        help="Fit with alpha")
    parser.add_argument(
        "--out", type=str, default="res.txt",
        help="Output file")
    parser.add_argument(
        "--outdir", type=str, default=".",
        help="Directory for output file")
    args = parser.parse_args()
   
    # Parse arguments =========================================================
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    resdir = args.resdir
    fixscale = args.fixscale
    chi2_min0 = args.chi2_min0
    # Parse arguments =========================================================


    # Read files and extract physical parameters and chi2 =====================
    # Note: All TPM results should be from the same conditions 
    #       (obs file, ephem file) with different parameters (A, TI, etc.).
    # Note: Input data should be homogeneous
    #       N_data = N_TI x N_HA (Haple angle, roughness) x N_A (albedo)
    # Calculate N_data, N_TI, N_HA, N_A, and N_obs
    A_list, TI_list, CA_list, CR_list = [], [], [], []
    scale_list, Dv_list, Ds_list, chi2_list = [], [], [], []
    Htheta_list = []
    N_obs = 0
    for idx_res, res in enumerate(args.res):
        #respath = os.path.join(resdir, res)
        respath = res
        with open (respath, "r") as f:
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
                    # Haple angle in deg (t.Bar in the output) 
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
                    # chi2 (not rduced yet)
                    chi2 = l_split[10]
                    chi2_list.append(float(chi2))
                elif (idx_res == 0) and ("ReadObs: Observation" in l):
                    # This selection is valid for both TI==0 and TI!=0
                    # Extract Ndata
                    l_split = l.split()
                    n_obs = int(l_split[4].split("=")[1])
                    N_obs += n_obs
                    
            print(f"    Physical parameters extracted.")

    print("")
    print(f"Number of observations is {N_obs}")
    print("")
    # Count nan

    print("Number of nan.")
    # A, TI, CA, CR are fixed parameters. (i.e., always wo/ nan)
    # N of nan is always the same maybe?
    nan_sc   = np.isnan(scale_list).sum()
    nan_Dv   = np.isnan(Dv_list).sum()
    nan_Ds   = np.isnan(Ds_list).sum()
    nan_chi2 = np.isnan(chi2_list).sum()
    print(f"    nan in sc   : {nan_sc}")
    print(f"    nan in Dv   : {nan_Dv}")
    print(f"    nan in Ds   : {nan_Ds}")
    print(f"    nan in chi2 : {nan_chi2}")
    # Read files and extract physical parameters and chi2 =====================
 

    N_all = len(A_list)
    TI_list = list(set(TI_list))
    TI_list = sorted(TI_list)
    N_TI = len(TI_list)

    A_list = list(set(A_list))
    A_list = sorted(A_list)
    N_A  = len(A_list)

    # This is the same with N_ca == N_cf (i.e., number of roughness parameters)
    Htheta_list = list(set(Htheta_list))
    N_Htheta  = len(Htheta_list)

    # This is also needed to identify filenames below
    CA_list = list(set(CA_list))
    N_CA  = len(CA_list)

    # Sort carefully......  !!!
    # Ascending order of Htheta and CA is not necessarily the same. 
    # (see, e.g., Hung+2022)
    CA_list = sorted(CA_list)
    Htheta_list = [crater2Htheta(ca, 999)for ca in CA_list]

    print(f"Input N_all      = {N_all}")
    print(f"      N_TI       = {N_TI}")
    print(f"      N_A        = {N_A}")
    print(f"      N_Htheta   = {N_Htheta}")
    print("")

    TI_thresh = args.TI_thresh
    TIrock_list = [x for x in TI_list if x >= TI_thresh]
    N_TIrock = len(TIrock_list)
    TIrego_list = [x for x in TI_list if x <= TI_thresh]
    N_TIrego = len(TIrego_list)
    print(f"Let's divide TI with TI_thresh = {TI_thresh}")
    print(f"      N_TI          = {N_TI}")
    print(f"      N_TIrock     = {N_TIrock}")
    print(f"      N_TIrego = {N_TIrego}")
    print(f"      Note: Not always N_TI = N_TIrock + N_TIrego")
    print("")

    # Regolith abundance 
    alpha_list = np.arange(args.amin, args.amax + args.astep, args.astep)
    print(f"List of regolith abundance: {alpha_list}")
    print("")
    # TODO: bond albedo, A, is always fixed?? For Eros, yes.
 

    # Calculate (chi2, alpha) for each (TI_rock, TI_rego, Htheta) =============
    N_comb = int(N_Htheta*N_TIrego*N_TIrock)
    idx_all = 1
    # Get all filenames
    resall = os.listdir(resdir)
    
    # Make DataFrame to register chi2
    column = ["Htheta", "TIrego", "TIrock", "alpha", "chi2"]
    df = pd.DataFrame(columns=column)


    # Loop for Htheta (or CA, i.e., roughness)
    for idx_Htheta, (Htheta, CA) in enumerate(zip(Htheta_list, CA_list)):
        # To extract two files
        # ex) tpmout_433_brute_A0.12_ti950_ca50_cr0.5.dat
        #  ->     *_ti{TIrego}_ca{CA}_*
        #     and *_ti{TIrock}_ca{CA}_*
        ca_str = f"_ca{int(CA)}_"

        # Loop for TI of regolith
        for idx_TIrego, TIrego in enumerate(TIrego_list):
            # Extract fluxes of regolith (frego)
            tirego_str = f"_ti{int(TIrego)}_"
            resrego = [x for x in resall if (tirego_str in x) and (ca_str in x)]
            frego = os.path.join(resdir, resrego[0])
            print(frego)
            df_rego = extract_flux(frego, fixscale)

            # Loop for TI of rock
            for idx_TIrock, TIrock in enumerate(TIrock_list):
                print(f"Htheta, CA, TIrego, TIrock = {Htheta}, {CA}, "
                      f"{TIrego}, {TIrock} ({idx_all}/{N_comb})")

                # Extract fluxes of rock (frock)
                tirock_str = f"_ti{int(TIrock)}_"
                resrock = [x for x in resall if (tirock_str in x) and (ca_str in x)]
                frock = os.path.join(resdir, resrock[0])
                df_rock = extract_flux(frock, fixscale)

                # Combine two dataframe and 
                #    return search best alpha and chi2 (if fitalpha == True)
                #    return all alpha and chi2         (if fitalpha == False)
                alpha_arr, chi2_arr = search_regolith_abundance(df_rego, df_rock, alpha_list, chi2_min0, args.fitalpha)
                
                # Save info.
                for a, c in zip(alpha_arr, chi2_arr):
                    print(f"  -> alpha, chi2 = {a:.2f}, {c:.2f}")
                    df.loc[len(df)] = [Htheta, TIrego, TIrock, a, c]

                # Update index
                idx_all += 1
    # Calculate (chi2, alpha) for each (TI_rock, TI_rego, Htheta) =============
     

    # Save
    df.to_csv(args.out, sep=" ", index=False)
