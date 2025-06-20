#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Search best regolith abundance etc. with NN using a brute-force method
with iterative processes.
The initial thermal inertia (TI0) is used to determine an initial cutoff (TI_thresh)
in the iterative processes.
Only the results when converged are saved as "blend_tpm_result.py".

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

from trem.emittance.common_emittance import (
    extract_flux, crater2Htheta, extract_bestparam)
from trem.emittance.common_dualcomponent import search_regolith_abundance
from trem.emittance.util_Cambioni2021 import calc_TIth


if __name__ == "__main__":
    parser = ap(
        description="Blend thermal fluxes introducing alpha, "
        "regolith abundance with iterative process")
    parser.add_argument(
        "res", type=str,
        help="Results of NN")
    parser.add_argument(
        "--TI0", type=float, default=150,
        help="Initial thermal inertia of rock to determine threshold of TI of regolith and rocks")
    parser.add_argument(
        "--TI_thresh", type=float, default=False,
        help="Thermal inertia cutoff")
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
        "--N_param", type=int, default=2,
        help="Number of parameters")
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
    fixscale = args.fixscale
    chi2_min0 = args.chi2_min0
    # Parse arguments =========================================================

    # Read files 
    df_NN = pd.read_csv(args.res, sep=" ")
    # Add dummy
    df_NN["scalefactor"] = 1
    Htheta_list = sorted(list(set(df_NN["Htheta"])))
    TI_list= sorted(list(set(df_NN["TI"])))
    df_temp = df_NN[
        (df_NN["Htheta"] == Htheta_list[0]) & (df_NN["TI"] == TI_list[0])
        ]
    N_Htheta = len(Htheta_list)
    N_TI = len(TI_list)

    # Number of data points
    N_data = len(df_temp)

    # Calculate degree of freedom
    N_param = args.N_param
    dof = N_data - N_param

    # Start iterative process =================================================
    print("Start iterative process to determine TI_threshold (a.k.a., Gamma_c)")
    print("")
    # Initial thermal inertia to determine threshold
    # (e.g., TI estimated with single-component TPM)
    TI_rock0 = args.TI0

    # Regolith abundance 
    alpha_list = np.arange(0, 1.0 + args.astep, args.astep)
    print(f"List of regolith abundance: {alpha_list}")
    print("")


    if args.TI_thresh:
        TI_thresh = args.TI_thresh

        # TI_thresh is given by hand
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

    else:
   
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
            # (See trem/emittance/util_Cambioni2021.py for detail)
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

 
            # Calculate (chi2, alpha) for each (TI_rock, TI_rego, Htheta) =====
            N_comb = int(N_Htheta*N_TIrego*N_TIrock)
            idx_all = 1
            
            # Make DataFrame to register chi2
            column = ["Htheta", "TIrego", "TIrock", "alpha", "chi2"]
            df = pd.DataFrame(columns=column)

            # Loop for Htheta (i.e., roughness)
            for idx_Htheta, Htheta, in enumerate(Htheta_list):
                print(f"  Htheta = {Htheta:.2f}")
                # Loop for TI of regolith
                for idx_TIrego, TIrego in enumerate(TIrego_list):
                    # Extract fluxes of regolith (frego)
                    df_rego = df_NN[(df_NN["Htheta"] == Htheta) & (df_NN["TI"] == TIrego)]
                    df_rego = df_rego.reset_index(drop=True)

                    # Loop for TI of rock
                    for idx_TIrock, TIrock in enumerate(TIrock_list):
                        df_rock = df_NN[(df_NN["Htheta"] == Htheta) & (df_NN["TI"] == TIrock)]
                        df_rock = df_rock.reset_index(drop=True)

                        # Combine two dataframe and return only alpha which gives
                        # the minimum chi2
                        # Note: results are already fit by alpha
                        alpha_arr, chi2_arr = search_regolith_abundance(
                            df_rego, df_rock, alpha_list, chi2_min0, True)
                        # Save info.
                        for a, c in zip(alpha_arr, chi2_arr):
                            # For test
                            #print(f"  -> alpha, chi2 = {a:.2f}, {c:.2f}")
                            df.loc[len(df)] = [Htheta, TIrego, TIrock, a, c]

                        # Update index
                        idx_all += 1
            # Calculate (chi2, alpha) for each (TI_rock, TI_rego, Htheta) =====
            
            # Determine the best fit parameters (TI_rock, TI_rego, Htheta)
            # Note: alpha is already fit)
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


    # Save the results when TI_th is converged.
    # (inherit TIrego_list and TIrock_list)
    # Make DataFrame to register chi2
    column = ["Htheta", "TIrego", "TIrock", "alpha", "chi2"]
    df = pd.DataFrame(columns=column)

    # Loop for Htheta (i.e., roughness)
    # (Maybe we can skip this 2nd calculation......, but I have no idea.)
    for idx_Htheta, Htheta, in enumerate(Htheta_list):
        print(f"  Htheta = {Htheta:.2f}")

        # Loop for TI of regolith
        for idx_TIrego, TIrego in enumerate(TIrego_list):
            # Extract fluxes of regolith (frego)
            df_rego = df_NN[(df_NN["Htheta"] == Htheta) & (df_NN["TI"] == TIrego)]
            df_rego = df_rego.reset_index(drop=True)

            # Loop for TI of rock
            for idx_TIrock, TIrock in enumerate(TIrock_list):
                # Extract fluxes of rock (frock)
                df_rock = df_NN[(df_NN["Htheta"] == Htheta) & (df_NN["TI"] == TIrock)]
                df_rock = df_rock.reset_index(drop=True)

                # Combine two dataframe and
                # return all alpha and chi2         (if fitalpha == False)
                alpha_arr, chi2_arr = search_regolith_abundance(
                    df_rego, df_rock, alpha_list, chi2_min0, False)
                # Save info.
                for a, c in zip(alpha_arr, chi2_arr):
                    #print(f"  -> alpha, chi2 = {a:.2f}, {c:.2f}")
                    df.loc[len(df)] = [Htheta, TIrego, TIrock, a, c]

    df.to_csv(args.out, sep=" ", index=False)
