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
import time
import numpy as np
import pandas as pd
from multiprocessing import Pool
from argparse import ArgumentParser as ap

from trem.common import elapsedtime
from trem.emittance.common_emittance import (
    extract_bestparam, extract_unique_epoch)
from trem.emittance.common_dualcomponent import (
    search_regolith_abundance, blend_flux_numpy)
from trem.emittance.util_Cambioni2021 import calc_TIth

def elapsedtime(t0, msg="Elapsed"):
    print(f"[TIME] {msg}: {time.time() - t0:.2f} s")

# --- Global cache ---
global_cache = None

def init_worker(shared_cache):
    global global_cache
    global_cache = shared_cache

def worker_Htheta_iter(args):
    (Htheta, TIrego_list, TIrock_list, alpha_list, chi2_min0) = args
    rows = []
    for TIrego in TIrego_list:
        df_rego = global_cache[(Htheta, TIrego)]
        for TIrock in TIrock_list:
            if TIrego > TIrock: continue
            df_rock = global_cache[(Htheta, TIrock)]
            alpha_arr, chi2_arr = search_regolith_abundance(df_rego, df_rock, alpha_list, chi2_min0, True)
            rows.extend([[Htheta, TIrego, TIrock, a, c] for a, c in zip(alpha_arr, chi2_arr)])
    return rows

def worker_Htheta_final(args):
    (Htheta, TIrego_list, TIrock_list, alpha_list, chi2_min0,
     fixscale, scale_all, scale_per_obs, out_file) = args

    rows_h = []
    for TIrego in TIrego_list:
        df_rego = global_cache[(Htheta, TIrego)]
        for TIrock in TIrock_list:
            if TIrego > TIrock: continue
            df_rock = global_cache[(Htheta, TIrock)]

            if fixscale:
                alpha_arr, chi2_arr = search_regolith_abundance(df_rego, df_rock, alpha_list, chi2_min0, False)
                rows_h.extend([[Htheta, TIrego, TIrock, a, c] for a, c in zip(alpha_arr, chi2_arr)])
            elif scale_all:
                for sf in np.arange(0.90, 1.11, 0.01):
                    d1, d2 = df_rego.copy(), df_rock.copy()
                    d1["scalefactor"], d2["scalefactor"] = sf, sf
                    alpha_arr, chi2_arr = search_regolith_abundance(d1, d2, alpha_list, chi2_min0, False)
                    rows_h.extend([[Htheta, TIrego, TIrock, a, c, sf] for a, c in zip(alpha_arr, chi2_arr)])
            elif scale_per_obs:
                t_unique, _ = extract_unique_epoch(df_rego, "jd")
                df_rego_w, df_rock_w = df_rego.copy(), df_rock.copy()
                for al in alpha_list:
                    sf_epoch_list = []
                    for epoch in t_unique:
                        r_ep, rk_ep = df_rego_w[df_rego_w["jd"]==epoch], df_rock_w[df_rock_w["jd"]==epoch]
                        best_sf, min_c = 1.0, float('inf')
                        for sf in np.arange(0.90, 1.10, 0.01):
                            f_b = blend_flux_numpy(r_ep["f_model"].to_numpy(), sf, rk_ep["f_model"].to_numpy(), sf, al)
                            chi2 = np.sum((r_ep["f_obs"].to_numpy()-f_b)**2 / r_ep["ferr_obs"].to_numpy()**2)
                            if chi2 < min_c: min_c, best_sf = chi2, sf
                        sf_epoch_list.append(best_sf)
                        df_rego_w.loc[df_rego_w["jd"]==epoch, "scalefactor"] = best_sf
                        df_rock_w.loc[df_rock_w["jd"]==epoch, "scalefactor"] = best_sf
                    f_b_f = blend_flux_numpy(df_rego_w["f_model"].to_numpy(), df_rego_w["scalefactor"].to_numpy(),
                                             df_rock_w["f_model"].to_numpy(), df_rock_w["scalefactor"].to_numpy(), al)
                    tot_c = np.sum((df_rego_w["f_obs"].to_numpy()-f_b_f)**2 / df_rego_w["ferr_obs"].to_numpy()**2)
                    rows_h.append([Htheta, TIrego, TIrock, al, tot_c] + sf_epoch_list)

    best_local, tmp_path = None, None
    if rows_h:
        if fixscale: cols = ["Htheta", "TIrego", "TIrock", "alpha", "chi2"]
        elif scale_all: cols = ["Htheta", "TIrego", "TIrock", "alpha", "chi2", "scalefactor"]
        else: cols = ["Htheta", "TIrego", "TIrock", "alpha", "chi2"] + [f"scalefactor{i+1}" for i in range(len(rows_h[0])-5)]

        df_chunk = pd.DataFrame(rows_h, columns=cols)
        tmp_path = f"{out_file}.tmp_{Htheta}.parquet"
        df_chunk.to_parquet(tmp_path, index=False)
        best_local = df_chunk.loc[df_chunk["chi2"].idxmin()].to_dict()
        del df_chunk
        rows_h.clear()
    return best_local, tmp_path

if __name__ == "__main__":
    parser = ap()
    parser.add_argument("res", type=str)
    parser.add_argument("--TI0", type=float, default=150)
    parser.add_argument("--bestparam", type=str, default=None)
    parser.add_argument("--TI_thresh", type=float, default=None)
    parser.add_argument("--obj", type=str, default="Eros")
    parser.add_argument("--T_typical", type=float, default=295.)
    parser.add_argument("--chi2_min0", type=float, default=200000)
    parser.add_argument("--phi", type=float, default=0.20)
    parser.add_argument("--astep", type=float, default=0.1)
    parser.add_argument("--fixscale", action="store_true")
    parser.add_argument("--scale_all", action="store_true")
    parser.add_argument("--scale_per_obs", action="store_true")
    parser.add_argument("--inbinary", action="store_true")
    parser.add_argument("--out", type=str, default="res.txt")
    parser.add_argument("--outbinary", action="store_true")
    parser.add_argument("--outsummary", type=str, default=None)
    args = parser.parse_args()

    t0 = time.time()

    # Read
    t_sub = time.time()
    df_NN = pd.read_parquet(args.res) if args.inbinary else pd.read_csv(args.res, sep=" ")
    df_NN["scalefactor"] = 1.0
    elapsedtime(t_sub, "Read CSV")

    # Extract
    t_sub = time.time()
    Htheta_list = sorted(df_NN["Htheta"].unique())
    TI_list = sorted(df_NN["TI"].unique())
    elapsedtime(t_sub, "Extract unique Htheta/TI")

    # Cache
    t_sub = time.time()
    df_cache = {k: v.reset_index(drop=True) for k, v in df_NN.groupby(["Htheta", "TI"])}
    del df_NN
    elapsedtime(t_sub, "Build cache")

    # --- Iterative process ---
    t_iter = time.time()
    # ここを TI_rock0 に統一
    TI_rock0 = pd.read_csv(args.bestparam, sep=" ")["TI"].iloc[0] if args.bestparam else args.TI0
    alpha_list = np.arange(0, 1.0 + args.astep, args.astep)

    if args.TI_thresh:
        TI_thresh = args.TI_thresh
        TIrk_l, TIrg_l = [x for x in TI_list if x >= TI_thresh], [x for x in TI_list if x <= TI_thresh]
    else:
        print("[INFO] Start iterative process to determine TI_thresh...")
        while True:
            _, TI_thresh = calc_TIth(TI_rock0, args.T_typical, args.obj, args.phi)
            TIrk_l, TIrg_l = [x for x in TI_list if x >= TI_thresh], [x for x in TI_list if x <= TI_thresh]

            print(f"  Current TI_rock0 = {TI_rock0:.2f}, Calculated TI_thresh = {TI_thresh:.2f}")
            print(f"  N_TIrock = {len(TIrk_l)}, N_TIrego = {len(TIrg_l)}")

            tasks = [(H, TIrg_l, TIrk_l, alpha_list, args.chi2_min0) for H in Htheta_list]
            t_p = time.time()
            with Pool(processes=4, initializer=init_worker, initargs=(df_cache,)) as pool:
                res = pool.map(worker_Htheta_iter, tasks)
            elapsedtime(t_p, "Parallel worker_Htheta_iter")

            df_it = pd.DataFrame([r for sub in res for r in sub], columns=["Htheta", "TIrego", "TIrock", "alpha", "chi2"])
            _, best_p = extract_bestparam(df_it, "chi2", ["TIrock"])
            TI_rock_best = best_p[0]
            dTI = abs(TI_rock_best - TI_rock0) / TI_rock0
            print(f"  Best TI_rock = {TI_rock_best:.2f}, dTI = {dTI:.5f}")

            if dTI < 1e-3: break
            TI_rock0 = TI_rock_best # 更新
    elapsedtime(t_iter, "Iterative TI_thresh")

    # --- Final calculation ---
    t_final = time.time()
    final_tasks = [(H, TIrg_l, TIrk_l, alpha_list, args.chi2_min0, args.fixscale, args.scale_all, args.scale_per_obs, args.out) for H in Htheta_list]

    t_p = time.time()
    with Pool(processes=4, initializer=init_worker, initargs=(df_cache,)) as pool:
        results = pool.map(worker_Htheta_final, final_tasks)
    elapsedtime(t_p, "Parallel worker_Htheta_final")

    # --- Result handling ---
    t_sub = time.time()
    valid_res = [r for r in results if r[1] is not None]
    if valid_res:
        best_final = pd.DataFrame([r[0] for r in valid_res]).sort_values("chi2").iloc[0]
        print(f"\n[SUMMARY] Best-fit parameters:")
        print(f"  chi2_min  = {best_final['chi2']:.2f}")
        print(f"  alpha     = {best_final['alpha']:.2f}")
        print(f"  Htheta    = {best_final['Htheta']:.2f}")
        print(f"  TIrego    = {best_final['TIrego']:.2f}")
        print(f"  TIrock    = {best_final['TIrock']:.2f}")

        df_final = pd.concat([pd.read_parquet(r[1]) for r in valid_res], ignore_index=True)
        elapsedtime(t_sub, "Build final DataFrame")

        t_save = time.time()
        is_binary = args.outbinary or args.out.endswith(".parquet")
        if is_binary:
            df_final.to_parquet(args.out, index=False)
        else:
            df_final.to_csv(args.out, sep=" ", index=False, float_format="%.2f")
        elapsedtime(t_save, "Save output file")

        for r in valid_res: os.remove(r[1]) # 掃除

    if args.outsummary:
        pd.DataFrame({"TI_thresh": [TI_thresh], "T_typical": [args.T_typical], "phi": [args.phi]}).to_csv(args.outsummary, sep=" ", index=False, float_format="%.2f")

    elapsedtime(t0, "Total")
    print()
