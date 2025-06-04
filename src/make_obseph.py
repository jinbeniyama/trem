#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Make obs file and ephem file for TPM.
Input argument, obs, should be dataframe(s), which consists of 
jd (light-time corrected!), wavelength, flux, fluxerr, code, cflag.

"""
from argparse import ArgumentParser as ap
import pandas as pd

from trem.common import make_ephemfile, make_obsfile


if __name__ == "__main__":
    parser = ap(description="Output obs file and ephem file for TPM.")
    parser.add_argument(
        "obj", type=str,
        help="object name")
    parser.add_argument(
        "obs", nargs="*", 
        help="jd, wavelength, flux, fluxerr, code, and cflag")
    parser.add_argument(
        "--rmnegativeflux", action="store_true", default=False,
        help="Remove negative fluxes and do not use them")
    parser.add_argument(
        "--ltcor", action="store_true", default=False,
        help="Do light-time correction in this code")
    parser.add_argument(
        "--out_obs", type=str, default="obs.txt",
        help="Obs file")
    parser.add_argument(
        "--out_eph", type=str, default="eph.txt",
        help="Ephemeris file")
    parser.add_argument(
        "--warmuptime_day", type=float, default=30.,
        help="Warming up time in day")
    args = parser.parse_args()

    # Merge all obs using jd
    # All inputs should have the same columns
    # jd, wavelength, flux, fluxerr, code, and cflag
    # Note: jd should be light-time corrected!
    if args.ltcor:
        print("  We do perform light-time correction in this script.")
        print("  Please check again that your jd is '''not light-time corrected!''' ")
    else:
        print("  We do not perform light-time correction in this script.")
        print("  Please check again that your jd is '''light-time corrected!''' ")
        print("  Vérifiez à nouveau que votre jd est corrigé en fonction de l'heure de la lumière!")
    
    # Read and merge observations
    df_list = []
    for x in args.obs:
        df = pd.read_csv(x, sep=" ")
        df_list.append(df)
    df_con = pd.concat(df_list)
    
    # Make ephem file
    make_ephemfile(args.obj, df_con, args.out_eph, args.warmuptime_day)

    # Make obs file w/ or wo/ light-time correction
    make_obsfile(args.obj, df_con, args.out_obs, args.ltcor, args.rmnegativeflux)
