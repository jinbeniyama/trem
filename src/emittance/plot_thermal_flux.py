#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot "only" thermal flux in outputs of TPM.
TODO: Plot optical flux as well.
"""
from argparse import ArgumentParser as ap

from tpmwrapper.common import extract_flux, plot_flux


if __name__ == "__main__":
    parser = ap(description="Plot thermal flux.")
    parser.add_argument(
        "res", nargs="*", 
        help="Results of rumtpm")
    parser.add_argument(
        "--fixscale", action="store_true",  default=False, 
        help="Fix scale parameter to 1")
    parser.add_argument(
        "--y1range", type=float, nargs=2, default=None, 
        help="Range of y1 axis (i.e., flux)")
    parser.add_argument(
        "--out", type=str, default="tflux.jpg",
        help="Output file")
    args = parser.parse_args()
   

    # Extract fluxes
    df_list = []
    for idx_res, res in enumerate(args.res):
        df = extract_flux(res, args.fixscale)
        df_list.append(df)

    # Plot residuals
    out = args.out
    plot_flux(df_list, args.fixscale, args.y1range, out)
    
