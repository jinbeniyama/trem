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
import sys
import time
from argparse import ArgumentParser as ap
import numpy as np
import pandas as pd

from trem.emittance.common_emittance import (
    extract_flux, crater2Htheta, extract_bestparam, introduce_var_scalefactor,
    extract_unique_epoch)
from trem.emittance.common_dualcomponent import search_regolith_abundance, blend_flux_numpy
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
        "--scale_all", action="store_true", default=False,
        help="Use global scale factor")
    parser.add_argument(
        "--scale_per_obs", action="store_true", default=False,
        help="Use scale factors per observation")
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
   
    t0 = time.time() 

    # Parse arguments =========================================================
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    fixscale = args.fixscale
    scale_all = args.scale_all
    scale_per_obs = args.scale_per_obs
    chi2_min0 = args.chi2_min0
    if fixscale & scale_per_obs:
        print("  Choose either fixscale or scale_per_obs.")
        print("  Exit")
        sys.exit()
    elif fixscale:
        print("  Scale factors are assumed to be 1.")
    elif scale_per_obs:
        print("  Scale factors are introduced per epoch. (only for spectroscopy)")
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
        print(f"Let's divide TI values into two with TI_thresh = {TI_thresh:.2f}")
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
                    # Extract fluxes of regolith
                    df_rego = df_NN[(df_NN["Htheta"] == Htheta) & (df_NN["TI"] == TIrego)]
                    df_rego = df_rego.reset_index(drop=True)

                    # Loop for TI of rock
                    for idx_TIrock, TIrock in enumerate(TIrock_list):
                        df_rock = df_NN[(df_NN["Htheta"] == Htheta) & (df_NN["TI"] == TIrock)]
                        df_rock = df_rock.reset_index(drop=True)

                        # Combine two dataframe and return only alpha which gives
                        # the minimum chi2
                        # Note: results are already fit by alpha
                        # Note: scale factors are assumed to be 1
                        # TODO: Should we introduce scale factors 
                        # as free parameters here as well?
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
    if fixscale:
        column = ["Htheta", "TIrego", "TIrock", "alpha", "chi2"]
    elif scale_all:
        # Note: This scale factor is applied for both spec. and phot.
        column = ["Htheta", "TIrego", "TIrock", "alpha", "chi2", "scalefactor"]
    elif scale_per_obs:
        pass


    # Loop for Htheta (i.e., roughness)
    # (Maybe we can skip this 2nd calculation......, but I have no idea.)
    
    rows_all = []
    for idx_Htheta, Htheta, in enumerate(Htheta_list):
        print(f"  Htheta = {Htheta:.2f}")

        # Loop for TI of regolith
        for idx_TIrego, TIrego in enumerate(TIrego_list):
            # Extract fluxes of regolith
            df_rego = df_NN[(df_NN["Htheta"] == Htheta) & (df_NN["TI"] == TIrego)].copy()
            df_rego = df_rego.reset_index(drop=True)

            # Loop for TI of rock
            for idx_TIrock, TIrock in enumerate(TIrock_list):
                # Extract fluxes of rock
                df_rock = df_NN[(df_NN["Htheta"] == Htheta) & (df_NN["TI"] == TIrock)].copy()
                df_rock = df_rock.reset_index(drop=True)

                # wo/ scale factors
                if fixscale:
                    # Combine two dataframe and return only alpha which gives
                    # the minimum chi2
                    # Note: results are already fit by alpha
                    # Note: scale factors are assumed to be 1
                    alpha_arr, chi2_arr = search_regolith_abundance(
                        df_rego, df_rock, alpha_list, chi2_min0, False)

                # w/ global scale factors
                elif scale_all:
                    # Combine two dataframe and return chi-squared values 

                    # TODO: As free parameters
                    sf0, sf1, sfstep = 0.90, 1.10, 0.01
                    sf_list = np.arange(sf0, sf1 + sfstep, sfstep)

                    key_t = "jd"
                    t_unique_list, _ = extract_unique_epoch(df_rego, key_t)
                    df_rego["scalefactor"] = df_rego["scalefactor"].astype(float)
                    df_rock["scalefactor"] = df_rock["scalefactor"].astype(float)

                    for sf in sf_list:
                        # Introduce scale factors for both spec. and phot.
                        df_rego.loc[:, "scalefactor"] = sf
                        df_rock.loc[:, "scalefactor"] = sf

                        sf_list1 = list(set(df_rego.scalefactor))
                        #print(f"  Unique scale factors: {sf_list1}")
                        alpha_arr, chi2_arr = search_regolith_abundance(
                            df_rego, df_rock, alpha_list, chi2_min0, False)
                        # Save info.
                        rows = [
                            [Htheta, TIrego, TIrock, a, c, sf]
                            for a, c in zip(alpha_arr, chi2_arr)
                        ]
                        rows_all.extend(rows)

                # w/ scale factors for spectra (not for photometry)
                elif scale_per_obs:
                    # Combine two dataframe and return chi-squared values 
                    # Note: Scale factors are introduced.
                    #       Results are already fit by the scale factors per obs.

                    # TODO: As free parameters
                    sf0, sf1, sfstep = 0.90, 1.10, 0.01
                    sf_list = np.arange(sf0, sf1, sfstep)

                    key_t = "jd"
                    t_unique_list, dfs_phot = extract_unique_epoch(df_rego, key_t)
                    df_rego["scalefactor"] = df_rego["scalefactor"].astype(float)
                    df_rock["scalefactor"] = df_rock["scalefactor"].astype(float)

                    # Test
                    alpha_list = [0]

                    # Search best scale parameters for each alpha
                    for al in alpha_list:

                        sf_epoch_list = []
                        for epoch in t_unique_list: 
                            df_rego_epoch = df_rego[df_rego["jd"] == epoch]
                            df_rock_epoch = df_rock[df_rock["jd"] == epoch]
                            
                            # Fit scale factor here
                            for idx_sf, sf in enumerate(sf_list):
                                # Introduce scale factors only for photometry
                                df_rego_epoch.loc[:, "scalefactor"] = sf
                                df_rock_epoch.loc[:, "scalefactor"] = sf

                                f1 = df_rego_epoch["f_model"].to_numpy()
                                f2 = df_rock_epoch["f_model"].to_numpy()
                                f_obs = df_rego_epoch["f_obs"].to_numpy()
                                ferr_obs = df_rego_epoch["ferr_obs"].to_numpy()

                                # Blended model flux for a combination of 
                                # (epoch, scale factor, alpha)
                                f_blend = blend_flux_numpy(f1, sf, f2, sf, al)

                                # Calculate chi2 of blended flux
                                diff = (f_obs - f_blend)**2 / ferr_obs**2
                                chi2 = np.sum(diff)
                                if idx_sf == 0:
                                    chi2_min_epoch_sf = chi2
                                    sf_epoch = sf
                                else:
                                    if chi2 < chi2_min_epoch_sf:
                                        chi2_min_epoch_sf = chi2
                                        sf_epoch = sf
                            #print(f"Best sf at {epoch} with alpha of {al}: {sf_epoch}")
                            # Update scale factors
                            df_rego.loc[df_rego["jd"]==epoch, "scalefactor"] = sf_epoch
                            df_rock.loc[df_rock["jd"]==epoch, "scalefactor"] = sf_epoch

                            sf_epoch_list.append(sf_epoch)

                        f1 = df_rego["f_model"].to_numpy()
                        sf_per_obs = df_rego["scalefactor"].to_numpy()
                        f2 = df_rock["f_model"].to_numpy()
                        f_obs = df_rock["f_obs"].to_numpy()
                        ferr_obs = df_rock["ferr_obs"].to_numpy()
                        f_blend = blend_flux_numpy(f1, sf_per_obs, f2, sf_per_obs, al)

                        # Calculate chi2 of blended flux
                        diff = (f_obs - f_blend)**2 / ferr_obs**2
                        chi2 = np.sum(diff)

                        # Save info.
                        # sf_epoch_list: best scale parameters for each epoch
                        rows = [[Htheta, TIrego, TIrock, al, chi2] + sf_epoch_list]
                        rows_all.extend(rows)
    if scale_per_obs:
        # Save all scale factors.
        column = ["Htheta", "TIrego", "TIrock", "alpha", "chi2"]
        for idx, epoch in enumerate(t_unique_list): 
            column.append(f"scalefactor{idx+1}")

    df = pd.DataFrame(rows_all, columns=column)

    df.to_csv(args.out, sep=" ", index=False, float_format="%.2f")
    t1 = time.time() 
    elapsed_time = t1 - t0
    print(f"  Elapsed time：{elapsed_time:.2f}")
