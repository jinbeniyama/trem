#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Do NEATM, not TPM, with various parameters, and plot the results.
Note: Input albedo of NEATM is "geometric albedo (pV)", not "Bond albedo (A)",
      which is input of TPM. See readme.pdf of the Marco's code.
"""
import os
from argparse import ArgumentParser as ap
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt  
import subprocess

from myplot import mycolor


if __name__ == "__main__":
    parser = ap(description="Do NEATM with various parameters.")
    parser.add_argument(
        "obj", type=str,
        help="object name")
    parser.add_argument(
        "inp", type=str,  
        help="Input file of fittm with nominal values")
    parser.add_argument(
        "--Herr", type=float, default=0.3,
        help="Uncertainty of H")
    parser.add_argument(
        "--Hstep", type=float, default=0.05,
        help="Step size of H")
    parser.add_argument(
        "--out", type=str, default="res.jpg",
        help="Output filename")
    args = parser.parse_args()
    

    Herr, Hstep = args.Herr, args.Hstep

    with open(args.inp, "r") as f:
        lines = f.readlines()
        # TODO: Update
        # This code works well for just one flux
        l0_s = lines[0].split()
        model = l0_s[0]
        H0    = float(l0_s[1])
        G     = l0_s[2]
        eps   = l0_s[3]
        # Not used.
        pv    = l0_s[4]
        eta   = l0_s[5]
        r     = l0_s[6]
        delta = l0_s[7]
        alpha = l0_s[8]
        # 
        l_flux = lines[1]
    
    # Do NEAM
    # geometric albedo
    pmin, pmax, pstep = 0.05, 0.70, 0.05
    pv_list = np.arange(pmin, pmax, pstep)
    # absolute magnitude
    H_list = np.arange(H0-Herr, H0+Herr, Hstep)
    inpdir = "inpfiles"
    
    # To save input
    pv0_list, H0_list = [], []
    Dkm1_list, pv1_list, eta1_list = [], [], []
    for pv in pv_list:
        for H in H_list:

            file_temp = f"{inpdir}/inp_{pv}_{H}.txt"
            with open(file_temp, "w") as file_temp:
                file_temp.write(f"{model} {H} {G} {eps} {eta} {pv} {r} {delta} {alpha}\n")
                file_temp.write(l_flux)
            file_temp = f"{inpdir}/inp_{pv}_{H}.txt"

            with open(file_temp, "r") as file_in:
                process = subprocess.Popen(["fittm"], stdin=file_in, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                stdout, stderr = process.communicate()
                print(stdout)
                
                # stdout is as below
                # l>    D(km)   sigmaD     pV    sigma_pV	 eta   sigma_eta   chi^2
                # o>    0.015   0.004	 0.249 0.097	 0.700 0.743	   0.0
                stdout = stdout.split("\n")
                res = stdout[1]
                res = res.split(" ")
                res = [x for x in res if x!=""]
                # Finally as below.
                # ['o>', '0.015', '0.006', '0.249', '0.138', '0.700', '0.805', '0.0']
                res = [x.replace("\t", "") for x in res]

                # Extract physical parameters
                Dkm1_list.append(float(res[1]))
                pv1_list.append(float(res[3]))
                eta1_list.append(float(res[5]))
                # Save input 
                pv0_list.append(pv)
                H0_list.append(H)

    df = pd.DataFrame(
        dict(H0=H0_list, pv0=pv0_list,
            Dkm1=Dkm1_list, pv1=pv1_list, eta1=eta1_list))

    # Save common info.
    df["obj"]    = args.obj
    df["Dm1"] = df["Dkm1"]*1e3
    

    # Plot results ============================================================
    fig = plt.figure(figsize=(16, 12))
    # Upper: inputs
    # Lower: outputs
    xwi, ywi = 0.22, 0.22
    colin, colout = "black", mycolor[0]

    # pv
    ax_u1 = fig.add_axes([0.10, 0.70, xwi, ywi])
    ax_u1.set_xlabel("input pv")
    ax_u1.set_ylabel("N")
    ax_u1.hist(
        df["pv0"], histtype="barstacked", bins=len(set(pv_list)), rwidth=0.6, color=colin)
    # H
    ax_u2 = fig.add_axes([0.40, 0.70, xwi, ywi])
    ax_u2.set_xlabel("input H")
    ax_u2.set_ylabel("N")
    ax_u2.hist(
        df["H0"], histtype="barstacked", bins=len(set(H_list)), rwidth=0.6, color=colin)
    # 
    #ax_u3 = fig.add_axes([0.70, 0.70, xwi, ywi])
    #ax_u3.set_xlabel("input H")
    #ax_u3.set_ylabel("N")
    #ax_u3.hist(df["H0"])

    # D
    nbin = 30
    ax_m1 = fig.add_axes([0.10, 0.40, xwi, ywi])
    ax_m1.set_xlabel("output D [m]")
    ax_m1.set_ylabel("N")
    ax_m1.hist(df["Dm1"], histtype="barstacked", bins=nbin, rwidth=0.6, color=colout)
    # pv
    ax_m2 = fig.add_axes([0.40, 0.40, xwi, ywi])
    ax_m2.set_xlabel("output pv")
    ax_m2.set_ylabel("N")
    ax_m2.hist(df["pv1"], histtype="barstacked", bins=nbin, rwidth=0.6, color=colout)
    # eta
    ax_m3 = fig.add_axes([0.70, 0.40, xwi, ywi])
    ax_m3.set_xlabel("output eta")
    ax_m3.set_ylabel("N")
    ax_m3.hist(df["eta1"], histtype="barstacked", bins=nbin, rwidth=0.6, color=colout)
    

    # D vs. pv
    ax_l1 = fig.add_axes([0.10, 0.1, xwi, ywi])
    ax_l1.set_xlabel("D [m]")
    ax_l1.set_ylabel("pv")
    ax_l1.scatter(df["Dm1"], df["pv1"], color=colout)
    # eta vs. pv
    ax_l2 = fig.add_axes([0.40, 0.1, xwi, ywi])
    ax_l2.set_xlabel("eta")
    ax_l2.set_ylabel("pv")
    ax_l2.scatter(df["eta1"], df["pv1"], color=colout)
    # D vs. eta
    ax_l3 = fig.add_axes([0.70, 0.1, xwi, ywi])
    ax_l3.set_xlabel("D [m]")
    ax_l3.set_ylabel("eta")
    ax_l3.scatter(df["Dm1"], df["eta1"], color=colout)
    fig.savefig(args.out)
    # Plot results ============================================================
