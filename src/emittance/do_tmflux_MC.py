#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Do tmflux with Monte-Carlo technique.

1. Sample (H, pV, eta) from "appropriate" distributions of them.
2. Calculate D with the equation:
   D = 1329./np.sqrt(pv) * 10**(-H/5)
3. Sample eta from normal distribution, the mean is determined 
   using empirical equation with phase angle:
    eta = 0.01*alpha +0.87.
4. Calculate flux F_8.7, F_10.7, and F_11.7 with NEATM.
5. N_mc sets of (H, pV, eta, D, F_8.7, F_10.7, F_11.7) are saved.
""" 
import numpy as np
import pandas as pd
from argparse import ArgumentParser as ap
import subprocess


def alpha2eta(alpha):
    """
    Convert phase angle to a typical eta.
    See Trilling+2016.

    Parameter
    ---------
    alpha : float
        phase angle
    
    Return
    ------
    eta : float
        beaming parameter
    """
    eta = 0.01*alpha +0.87
    return eta


def pvWright2016(N, seed=0):
    """
    Sample N pv from double peak pv distribution in Wright+2016.
    Note that the distribution is made using H < 22 NEAs.

    Parameters
    ----------
    N : int
        number of samples
    seed : int, optional
        random seed

    Return
    ------
    pv_list : array-like
        list of albedo with a length of N
    """
    np.random.seed(seed)
    pv_list = []
    for n in range(N):
        # Two random numbersa
        x = np.random.rand()
        y = np.random.rand() 
        # Dark fraction, peak of bright pv, peak of dark pv
        fD, b, d = 0.253, 0.168, 0.030
        if x < fD:
            t = d
        else:
            t = b
        pv = t*np.sqrt(-2*np.log(1-y))
        pv_list.append(pv)
    return pv_list


if __name__ == "__main__":
    parser = ap(description="Do tmflux with Monte Carlo method.")
    parser.add_argument(
        "--obj", type=str, default="1998 KY26",
        help="Object")
    parser.add_argument(
        "--etadist", type=str, default="normal",
        help="Eta distribution (normal or uniform)")
    parser.add_argument(
        "--Nmc", type=int, default=1000,
        help="Number of trials")
    parser.add_argument(
        "--out", type=str, default="NEATM_MCres.txt",
        help="Output file name")
    args = parser.parse_args()

    N_mc = args.Nmc 
    obj = args.obj

    # 2024-05-26
    if obj == "1998 KY26":
        r, delta, alpha = 1.04, 0.036, 41.4
        H, Herr = 26.3, 0.3
        # ?
        G, Gerr = 0.15, 0.10
    elif obj == "2009 BD":
    # 2013-10-13, 14
        r, delta, alpha = 1.075, 1.455, 56.5
        # Micheli+2012, which is also used in Mommert2016
        H, Herr = 28.43, 0.12
        G, Gerr = 0.18, 0.13
    else:
        r, delta, alpha = args.r, args.delta, args.alpha
        H, Herr = args.H, args.Herr
        G, Gerr = args.G, args.Gerr

    print(f"Input parameters:")
    print(f"  r, delta, alpha = ({r} au, {delta} au, {alpha} deg)")
    print(f"               H  = {H}+-{Herr}")
    print(f"               G  = {G}+-{Gerr}")

    # From Wright+2016
    pv_list = pvWright2016(N_mc)
    H_list = np.random.normal(H, Herr, N_mc)
    G_list = np.random.normal(G, Gerr, N_mc)

    # Determin eta from phase angle
    if args.etadist == "normal":
        # From Trilling+2016
        eta = alpha2eta(alpha)
        w_eta = 0.5
        eta_list = np.random.normal(eta, w_eta, N_mc)
    elif args.etadist == "uniform":
        # Uniform eta
        eta_list =np.random.uniform(0, 3, N_mc)

    
    D_list = []
    f4_5_list = []
    f8_7_list = []
    f10_7_list = []
    f11_7_list = []
    for idx in range(N_mc):
        for percent in np.arange(10, 101, 10):
            if (idx+1) == int(N_mc*percent/100):
                print(f"  {percent} percent finished {idx+1}/{N_mc}")
        D = 1329./np.sqrt(pv_list[idx]) * 10**(-H_list[idx]/5)
        #cmd = "echo "+str(D)+" "+str(pv_list[idx])+" " + str(eta_list[idx]) + \
        #f" {r} {delta} {alpha} | tmflux -l 10.7 11 1"    
        
        # TODO: emissivity?
        arg = f"{D} {pv_list[idx]} {eta_list[idx]} {r} {delta} {alpha} | tmflux -r -G {G_list[idx]}"
        
        # Extract B10.7
        pp = subprocess.Popen(
            ["tmflux", "-l", "4.5",  "11.9",  "0.2"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout, stderr = pp.communicate(arg)
        # This is like
        #   ['8.700', '9.555779e-04', '9.700', '1.178966e-03', '10.700', '1.363904e-03', '11.700', '1.507815e-03']
        f = stdout.split()
        # Extract 4.5 and 8.7, 10.7, and 11.7
        f4_5  = float(f[1])
        f8_7  = float(f[43])
        f10_7 = float(f[63])
        f11_7 = float(f[73])
        
        # In Jy
        f4_5_list.append(f4_5)
        f8_7_list.append(f8_7)
        f10_7_list.append(f10_7)
        f11_7_list.append(f11_7)
        D_list.append(D)

    
    df = pd.DataFrame(
        dict(H=H_list, G=G_list, pv=pv_list, D=D_list, eta=eta_list,  
        flux4=f4_5_list, flux8=f8_7_list, flux10=f10_7_list, flux11=f11_7_list))
    df.to_csv(args.out, sep=" ")
