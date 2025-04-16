#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Search best regolith abundance etc. using dual-component & brute-force method
with iterative processes.
The initial thermal inertia (TI0) is used to determine an initial cutoff (TI_thresh)
in the iterative processes.
Only the results when converged are saved as "blend_tpm_result.py".

TODO:
Scale parameter per spectrum.

Example
-------
# Calculate best chi2 values with variety of parameters.
> blend_tpm_result_iter.py tpmout* 

TODO
----
bond albedo, A, is always fixed?? For Eros, yes.
"""
import os 
from argparse import ArgumentParser as ap
import numpy as np
import pandas as pd

from tpmwrapper.common import (
    extract_flux, search_regolith_abundance, crater2Htheta, extract_bestparam)
from tpmwrapper.util_Cambioni2021 import calc_TIth


if __name__ == "__main__":
    parser = ap(
        description="Blend thermal fluxes introducing alpha, "
        "regolith abundance with iterative process")
    parser.add_argument(
        "res", nargs="*", 
        help="Results of rumtpm")
    parser.add_argument(
        "--resdir", type=str, default=".",
        help="Directory in which results are saved")
    parser.add_argument(
        "--TI0", type=float, default=150,
        help="Initial thermal inertia of rock to determine threshold of TI of regolith and rocks")
    parser.add_argument(
        "--obj", type=str, default="Eros",
        help="Object to refer to physical parameters")
    parser.add_argument(
        "--T_typical", type=float, default=295.,
        help="Typical temperature in K")
    parser.add_argument(
        "--chi2_min0", type=float, default=200000,
        help="Initial minimum chi2 used to find the best alpha")
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
    Htheta_list = [crater2Htheta(ca, 999) for ca in CA_list]

    print(f"Input N_all      = {N_all}")
    print(f"      N_TI       = {N_TI}")
    print(f"      N_A        = {N_A}")
    print(f"      N_Htheta   = {N_Htheta}")
    print("")




    # Start iterative process =================================================
    print("Start iterative process to determine TI_threshold (a.k.a., Gamma_c)")
    print("")
    # Initial thermal inertia to determine threshold
    # (e.g., TI estimated with single-component TPM)
    TI_rock0 = args.TI0

    # Regolith abundance 
    alpha_list = np.arange(0, 1.0+args.astep, args.astep)
    print(f"List of regolith abundance: {alpha_list}")
    print("")
   
    # Iterate until TI_th converges (dTI < dTI_goal)
    dTI  = 1e5
    dTI_goal = 1e-3
    while True:
        # TODO:
        # How to determine a typical temperature (T_typical)?
        # To determine c_p(T).
        # TI prop c_p**0.5, so the dependence is usually weak.
        # But it is better to consider later.
        T_typical = args.T_typical
        # Determine a TI_threshold with TI_rock
        # (See util_Cambioni2021.py for detail)
        _, TI_thresh = calc_TIth(TI_rock0, T_typical, args.obj)

        # Make lists of TIrock and TIrego with a given TI_thresh
        TIrock_list = [x for x in TI_list if x >= TI_thresh]
        N_TIrock = len(TIrock_list)
        TIrego_list = [x for x in TI_list if x <= TI_thresh]
        N_TIrego = len(TIrego_list)
        print(f"Let's divide TI with TI_thresh = {TI_thresh:.2f}")
        print(f"      N_TI          = {N_TI}")
        print(f"      N_TIrock      = {N_TIrock}")
        print(f"      N_TIrego      = {N_TIrego}")
        print(f"      Note: Not always N_TI = N_TIrock + N_TIrego")
        print("")

 
        # Calculate (chi2, alpha) for each (TI_rock, TI_rego, Htheta) =========
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
                df_rego = extract_flux(frego, fixscale)

                # Loop for TI of rock
                for idx_TIrock, TIrock in enumerate(TIrock_list):
                    #print(f"Htheta, CA, TIrego, TIrock = {Htheta}, {CA}, "
                    #      f"{TIrego}, {TIrock} ({idx_all}/{N_comb})")

                    # Extract fluxes of rock (frock)
                    tirock_str = f"_ti{int(TIrock)}_"
                    resrock = [x for x in resall if (tirock_str in x) and (ca_str in x)]
                    frock = os.path.join(resdir, resrock[0])
                    df_rock = extract_flux(frock, fixscale)

                    # Combine two dataframe and 
                    #    return search best alpha and chi2 (with fitalpha == True)
                    # Results are already fit by alpha
                    alpha_arr, chi2_arr = search_regolith_abundance(
                        df_rego, df_rock, alpha_list, chi2_min0, True)
                    # Save info.
                    for a, c in zip(alpha_arr, chi2_arr):
                        # For test
                        #print(f"  -> alpha, chi2 = {a:.2f}, {c:.2f}")
                        df.loc[len(df)] = [Htheta, TIrego, TIrock, a, c]

                    # Update index
                    idx_all += 1
        # Calculate (chi2, alpha) for each (TI_rock, TI_rego, Htheta) =========
        

        # Determine the best fit parameters (note: alpha is already fit)
        # TI_rock, TI_rego, Htheta
        key_chi2 = "chi2"
        # Use only TIrock to check the convergence
        params = ["TIrock"]
        chi2_min, best_params = extract_bestparam(df, key_chi2, params)
        TI_rock_best = best_params[0]

        # Check convergence
        dTI = abs(TI_rock_best - TI_rock0)/TI_rock0
        print(f"      -> TI_rock0       = {TI_rock0:.2f}")
        print(f"         TI_rock_best   = {TI_rock_best:.2f}")
        print(f"         dTI            = {dTI:.2f}")

        # Finish the iterative process
        if dTI < dTI_goal:
            print(f"         Converged.")
            print("")
            break
        # Otherwise update initial thermal inertia of rock (TI_rock)
        else:
            TI_rock0 = TI_rock_best
            print(f"         Not converged yet.")
            print("")

    
    # Save the results when TI_th is converged
    # Make DataFrame to register chi2
    column = ["Htheta", "TIrego", "TIrock", "alpha", "chi2"]
    df = pd.DataFrame(columns=column)

    # Loop for Htheta (or CA, i.e., roughness)
    # (Maybe we can skip this 2nd calculation......, but I have no idea.)
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
            df_rego = extract_flux(frego, fixscale)

            # Loop for TI of rock
            for idx_TIrock, TIrock in enumerate(TIrock_list):
                #print(f"Htheta, CA, TIrego, TIrock = {Htheta}, {CA}, "
                #      f"{TIrego}, {TIrock} ({idx_all}/{N_comb})")

                # Extract fluxes of rock (frock)
                tirock_str = f"_ti{int(TIrock)}_"
                resrock = [x for x in resall if (tirock_str in x) and (ca_str in x)]
                frock = os.path.join(resdir, resrock[0])
                df_rock = extract_flux(frock, fixscale)

                # Combine two dataframe and
                # return all alpha and chi2         (if fitalpha == False)
                alpha_arr, chi2_arr = search_regolith_abundance(
                    df_rego, df_rock, alpha_list, chi2_min0, False)
                # Save info.
                for a, c in zip(alpha_arr, chi2_arr):
                    #print(f"  -> alpha, chi2 = {a:.2f}, {c:.2f}")
                    df.loc[len(df)] = [Htheta, TIrego, TIrock, a, c]

    df.to_csv(args.out, sep=" ", index=False)
