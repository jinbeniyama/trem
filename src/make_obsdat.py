#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Make obs.dat before making obs file and ephem file for TPM.
The format is as follows.
``
jd wavelength flux fluxerr code cflag memo
2452545.9060276826 8.12 6.290097664831847 0.07037871513098569 568 999 UKIRT2002
``
"""
from argparse import ArgumentParser as ap
import pandas as pd


if __name__ == "__main__":
    parser = ap(description="Output obs.dat.")
    parser.add_argument(
        "jd0", type=float, 
        help="starting time of observations")
    parser.add_argument(
        "jd1", type=float,
        help="Ending time of observations")
    parser.add_argument(
        "jdstep", type=float,
        help="Time step in jd")
    parser.add_argument(
        "--code", type=str, default="500",
        help="MPC observatory code")
    parser.add_argument(
        "--wavelength", nargs="*", default=[8.7],
        help="Wavelength of interest")
    parser.add_argument(
        "--out", type=str, default="obs.dat",
        help="Obs.dat file")
    args = parser.parse_args()
    


    jd0, jd1, jdstep = args.jd0, args.jd1, args.jdstep
    jd = jd0
    w_list = args.wavelength
    code = args.code
    
    with open(args.out, "w") as f:
        head = "jd wavelength flux fluxerr code cflag memo\n"
        f.write(head)
        while jd < jd1:
            for w in w_list:
                line = f"{jd} {w} {999} {999} {code} 0 NoMemo\n"
                f.write(line)
            jd += jdstep
