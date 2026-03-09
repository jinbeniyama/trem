#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Search best regolith abundance etc. with NN using a brute-force method
with iterative processes.

Optimized version:
- Parallelized over Htheta
- Cached dataframe filtering
"""

import os
import sys
import time
from argparse import ArgumentParser as ap
from multiprocessing import Pool, cpu_count

import numpy as np
import pandas as pd

from trem.common import elapsedtime
from trem.emittance.common_emittance import (
    extract_bestparam, extract_unique_epoch)
from trem.emittance.common_dualcomponent import (
    search_regolith_abundance, blend_flux_numpy)
from trem.emittance.util_Cambioni2021 import calc_TIth


# ---------------------------------------------------------
# Cache dataframe slices
# ---------------------------------------------------------

def build_cache(df_NN, Htheta_list, TI_list):

    cache = {}

    for Htheta in Htheta_list:
        for TI in TI_list:

            df_tmp = df_NN[
                (df_NN["Htheta"] == Htheta) &
                (df_NN["TI"] == TI)
            ].reset_index(drop=True)

            cache[(Htheta, TI)] = df_tmp

    return cache


# ---------------------------------------------------------
# Worker for multiprocessing
# ---------------------------------------------------------

def worker_Htheta_iter(args):

    (
        Htheta,
        TIrego_list,
        TIrock_list,
        alpha_list,
        chi2_min0,
        df_cache
    ) = args

    rows = []

    for TIrego in TIrego_list:

        df_rego = df_cache[(Htheta, TIrego)]

        for TIrock in TIrock_list:

            df_rock = df_cache[(Htheta, TIrock)]

            alpha_arr, chi2_arr = search_regolith_abundance(
                df_rego,
                df_rock,
                alpha_list,
                chi2_min0,
                True
            )

            for a, c in zip(alpha_arr, chi2_arr):
                rows.append([Htheta, TIrego, TIrock, a, c])

    return rows


# Inside worker function for final calculation
def worker_Htheta_final(args):

    (
        Htheta,
        TIrego_list,
        TIrock_list,
        alpha_list,
        chi2_min0,
        df_cache,
        fixscale,
        scale_all,
        scale_per_obs
    ) = args

    rows_all = []

    for TIrego in TIrego_list:

        df_rego = df_cache[(Htheta, TIrego)].copy()

        for TIrock in TIrock_list:

            # Skip meaningless combination
            if TIrego > TIrock:
                continue

            df_rock = df_cache[(Htheta, TIrock)].copy()

            if fixscale:

                alpha_arr, chi2_arr = search_regolith_abundance(
                    df_rego,
                    df_rock,
                    alpha_list,
                    chi2_min0,
                    False
                )

                rows = [
                    [Htheta, TIrego, TIrock, a, c]
                    for a, c in zip(alpha_arr, chi2_arr)
                ]

                rows_all.extend(rows)

            elif scale_all:

                sf0, sf1, sfstep = 0.90, 1.10, 0.01
                sf_list = np.arange(sf0, sf1 + sfstep, sfstep)

                for sf in sf_list:

                    df_rego.loc[:, "scalefactor"] = sf
                    df_rock.loc[:, "scalefactor"] = sf

                    alpha_arr, chi2_arr = search_regolith_abundance(
                        df_rego,
                        df_rock,
                        alpha_list,
                        chi2_min0,
                        False
                    )

                    rows = [
                        [Htheta, TIrego, TIrock, a, c, sf]
                        for a, c in zip(alpha_arr, chi2_arr)
                    ]

                    rows_all.extend(rows)

            elif scale_per_obs:

                sf0, sf1, sfstep = 0.90, 1.10, 0.01
                sf_list = np.arange(sf0, sf1, sfstep)

                key_t = "jd"

                t_unique_list, _ = extract_unique_epoch(df_rego, key_t)

                df_rego["scalefactor"] = df_rego["scalefactor"].astype(float)
                df_rock["scalefactor"] = df_rock["scalefactor"].astype(float)

                for al in alpha_list:

                    sf_epoch_list = []

                    for epoch in t_unique_list:

                        df_rego_epoch = df_rego[df_rego["jd"] == epoch]
                        df_rock_epoch = df_rock[df_rock["jd"] == epoch]

                        for idx_sf, sf in enumerate(sf_list):

                            df_rego_epoch.loc[:, "scalefactor"] = sf
                            df_rock_epoch.loc[:, "scalefactor"] = sf

                            f1 = df_rego_epoch["f_model"].to_numpy()
                            f2 = df_rock_epoch["f_model"].to_numpy()
                            f_obs = df_rego_epoch["f_obs"].to_numpy()
                            ferr_obs = df_rego_epoch["ferr_obs"].to_numpy()

                            f_blend = blend_flux_numpy(f1, sf, f2, sf, al)

                            diff = (f_obs - f_blend) ** 2 / ferr_obs ** 2
                            chi2 = np.sum(diff)

                            if idx_sf == 0:

                                chi2_min_epoch_sf = chi2
                                sf_epoch = sf

                            else:

                                if chi2 < chi2_min_epoch_sf:

                                    chi2_min_epoch_sf = chi2
                                    sf_epoch = sf

                        df_rego.loc[
                            df_rego["jd"] == epoch,
                            "scalefactor"
                        ] = sf_epoch

                        df_rock.loc[
                            df_rock["jd"] == epoch,
                            "scalefactor"
                        ] = sf_epoch

                        sf_epoch_list.append(sf_epoch)

                    f1 = df_rego["f_model"].to_numpy()
                    sf_per_obs = df_rego["scalefactor"].to_numpy()
                    f2 = df_rock["f_model"].to_numpy()
                    f_obs = df_rock["f_obs"].to_numpy()
                    ferr_obs = df_rock["ferr_obs"].to_numpy()

                    f_blend = blend_flux_numpy(
                        f1,
                        sf_per_obs,
                        f2,
                        sf_per_obs,
                        al
                    )

                    diff = (f_obs - f_blend) ** 2 / ferr_obs ** 2
                    chi2 = np.sum(diff)

                    rows = [[Htheta, TIrego, TIrock, al, chi2] + sf_epoch_list]

                    rows_all.extend(rows)

    return rows_all

# ---------------------------------------------------------
# Main
# ---------------------------------------------------------

if __name__ == "__main__":
    parser = ap()
    parser.add_argument(
        "res", type=str)
    parser.add_argument(
        "--TI0", type=float, default=150)
    parser.add_argument(
        "--TI_thresh", type=float, default=False)
    parser.add_argument(
        "--obj", type=str, default="Eros")
    parser.add_argument(
        "--T_typical", type=float, default=295.)
    parser.add_argument(
        "--chi2_min0", type=float, default=200000)
    parser.add_argument(
        "--phi", type=float, default=0.20)
    parser.add_argument(
        "--astep", type=float, default=0.1)
    parser.add_argument(
        "--fixscale", action="store_true", default=False)
    parser.add_argument(
        "--scale_all", action="store_true", default=False)
    parser.add_argument(
        "--scale_per_obs", action="store_true", default=False)
    parser.add_argument(
        "--out", type=str, default="res.txt")
    args = parser.parse_args()

    t0 = time.time()

    df_NN = pd.read_csv(args.res, sep=" ")

    df_NN["scalefactor"] = 1

    Htheta_list = sorted(list(set(df_NN["Htheta"])))
    TI_list = sorted(list(set(df_NN["TI"])))

    df_cache = build_cache(df_NN, Htheta_list, TI_list)

    alpha_list = np.arange(0, 1.0 + args.astep, args.astep)

    # After parsing arguments
    TI_rock0 = args.TI0
    dTI_goal = 1e-3

    # If TI_thresh is given, skip iteration
    if args.TI_thresh:
        TI_thresh = args.TI_thresh
        print(f"[INFO] TI_thresh manually set: {TI_thresh:.2f}")

        TIrock_list = [x for x in TI_list if x >= TI_thresh]
        TIrego_list = [x for x in TI_list if x <= TI_thresh]

    else:
        # Iterative determination of TI_thresh
        print("[INFO] Start iterative process to determine TI_thresh...")
        while True:
            _, TI_thresh = calc_TIth(
                TI_rock0,
                args.T_typical,
                args.obj,
                args.phi
            )

            TIrock_list = [x for x in TI_list if x >= TI_thresh]
            TIrego_list = [x for x in TI_list if x <= TI_thresh]

            print(f"  Current TI_rock0 = {TI_rock0:.2f}")
            print(f"  Calculated TI_thresh = {TI_thresh:.2f}")
            print(f"  N_TIrock = {len(TIrock_list)}, N_TIrego = {len(TIrego_list)}")

            tasks = [
                (
                    Htheta,
                    TIrego_list,
                    TIrock_list,
                    alpha_list,
                    args.chi2_min0,
                    df_cache
                )
                for Htheta in Htheta_list
            ]

            with Pool(cpu_count()) as pool:
                results = pool.map(worker_Htheta_iter, tasks)

            rows = [r for sub in results for r in sub]

            df = pd.DataFrame(
                rows,
                columns=["Htheta", "TIrego", "TIrock", "alpha", "chi2"]
            )

            chi2_min, best_params = extract_bestparam(
                df,
                "chi2",
                ["TIrock"]
            )

            TI_rock_best = best_params[0]
            dTI = abs(TI_rock_best - TI_rock0) / TI_rock0

            print(f"  Best TI_rock = {TI_rock_best:.2f}, dTI = {dTI:.5f}")

            if dTI < dTI_goal:
                print(f"[INFO] Converged: TI_rock_best = {TI_rock_best:.2f}, TI_thresh = {TI_thresh:.2f}\n")
                break
            else:
                TI_rock0 = TI_rock_best
                print("  Not yet converged, updating TI_rock0...\n")

    
    tasks = [

        (
            Htheta,
            TIrego_list,
            TIrock_list,
            alpha_list,
            args.chi2_min0,
            df_cache,
            args.fixscale,
            args.scale_all,
            args.scale_per_obs
        )

        for Htheta in Htheta_list

    ]

    with Pool(cpu_count()) as pool:

        results = pool.map(worker_Htheta_final, tasks)

    rows_all = [r for sub in results for r in sub]
    
    # ---------------------------------------------------------
    # After final calculation and before saving CSV
    # ---------------------------------------------------------
    
    if args.fixscale or args.scale_all or args.scale_per_obs:
        # Convert results to DataFrame
        if args.fixscale:
            column = ["Htheta", "TIrego", "TIrock", "alpha", "chi2"]
        elif args.scale_all:
            column = ["Htheta", "TIrego", "TIrock", "alpha", "chi2", "scalefactor"]
        elif args.scale_per_obs:
            # For per-epoch scale factors, keep columns for each epoch
            column = ["Htheta", "TIrego", "TIrock", "alpha", "chi2"]
            for idx, epoch in enumerate(range(len(rows_all[0]) - 5)):  # adjust
                column.append(f"scalefactor{idx+1}")
    
    df = pd.DataFrame(rows_all, columns=column)
    
    # -----------------------------
    # Print summary of best chi2
    # -----------------------------
    best_idx = df["chi2"].idxmin()
    best_row = df.loc[best_idx]
    
    print("\n[SUMMARY] Best-fit parameters:")
    print(f"  chi2_min  = {best_row['chi2']:.2f}")
    print(f"  alpha     = {best_row['alpha']:.2f}")
    print(f"  Htheta    = {best_row['Htheta']:.2f}")
    print(f"  TIrego    = {best_row['TIrego']:.2f}")
    print(f"  TIrock    = {best_row['TIrock']:.2f}")
    
    if args.scale_all:
        print(f"  scalefactor = {best_row['scalefactor']:.3f}")
    elif args.scale_per_obs:
        sf_cols = [c for c in df.columns if "scalefactor" in c]
        sf_list = [best_row[c] for c in sf_cols]
        print(f"  scale factors per epoch = {sf_list}")
    
    # Save CSV
    df.to_csv(args.out, sep=" ", index=False, float_format="%.2f")
    
    elapsedtime(t0)
