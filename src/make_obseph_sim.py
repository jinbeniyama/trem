#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Make obs file and ephem file for TPM from simulation.

TODO: Add 1 at last?
"""
import os 
from argparse import ArgumentParser as ap


if __name__ == "__main__":
    parser = ap(description="Make obs and eph file for tpm.")
    parser.add_argument(
        "--obs", type=str, default="obs.txt",
        help="Output file")
    parser.add_argument(
        "--eph", type=str, default="eph.txt",
        help="Output file")
    parser.add_argument(
        "--rotP_s", type=float, default=60.0,
        help="Rotation period in s")
    parser.add_argument(
        "--t_sample", type=float, default=1.0,
        help="Sampling timescale in s")
    parser.add_argument(
        "--alpha", type=int, default=45,
        help="Phase angle in deg")
    parser.add_argument(
        "--N_rot", type=int, default=2,
        help="Number of rotations")
    parser.add_argument(
        "--outdir", type=str, default=".",
        help="Directory for output file")
    args = parser.parse_args()
   
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
     
    # Parse arguments
    N_rot = args.N_rot
    rotP_s = args.rotP_s
    rotP_hr = rotP_s/3600.
    rotP_d  = rotP_hr/24.
    alpha  = args.alpha

    # Related to samplling 
    dt_s = args.t_sample
    dt_d = dt_s/3600./24.


    # UPDATE ==================================================================
    # Very simple geometry
    
    # x_r, y_r, z_r
    # the position of the body in the heliocentric ecliptic reference frame 
    # as SEEN AT THE ASTEROID. (i.e., Just the location)
    # From asteroid to Sun (Ast -> Sun)
    
    # x_o, y_o, z_o
    # From asteroid to Earth
    # Sun --- Earth - Asteroid (definition)
    # From asteroid to Observer (Earth) (Ast -> Observer)


    # 1) alpha = 0 deg
    #
    # Sun --- (1 au) --- Earth --- (0.01 au) --- Ast 
    if alpha == 0:
        x_r, y_r, z_r = -1.01, 0.0, 0.0
        x_o, y_o, z_o = -0.01, 0.0, 0.0

    # 2) alpha ~ 30 deg
    #                            - - Ast
    # Sun --- (1 au) --- Earth - 
    elif alpha == 45:
        x_r, y_r, z_r = -1.01, -0.01, 0.0
        x_o, y_o, z_o = -0.01, -0.01, 0.0
    
    # 2) alpha ~ 60 deg
    #                              _ Ast
    #                            -  
    # Sun --- (1 au) --- Earth - 
    elif alpha == 60:
        x_r, y_r, z_r = -1.01,  -0.02, 0.0
        x_o, y_o, z_o = -0.01,  -0.02, 0.0
    else:
        assert False, "Not implemented."

    # UPDATE ==================================================================


    # Ephemeris file ==========================================================
    # JD r_X r_Y r_Z
    # the position of the body in the heliocentric ecliptic reference frame 
    # as SEEN AT THE ASTEROID. (i.e., Just the location)
    # Ast -> Sun vector

    # Consider Margin
    t0 = 30*rotP_d
    t = t0
    # N_rot rotation
    t1 = t0 + N_rot*rotP_d
   
    out = os.path.join(args.outdir, args.eph)
    Nsampling = 0
    with open(out, "w") as f:
        f.write(f"0 {x_r} {x_r} {x_r}\n")
        while t < t1:
            print(t)
            f.write(f"{t} {x_r} {x_r} {x_r}\n")
            t += dt_d
            Nsampling += 1
    # Ephemeris file ==========================================================


    # Observation file ========================================================
    # --------------
    # Nobs
    # 
    # JD Ndata
    # r_X r_Y r_Z (as for ephemeris file)
    # x y z
    # heliocentric X,Y,Z components of the vector from the asteroid to the Earth 
    # as SEEN AT THE ASTEROID. (i.e., Just the location)


    # Consider Margin
    t = t0
    # N_rot rotation
    t1 = t0 + N_rot*rotP_d
 
    # Nobs = Ndata*Nsampling
    # Only 10 micron
    Ndata = 1
    Nobs = Ndata*Nsampling

    out = os.path.join(args.outdir, args.obs)
    with open(out, "w") as f:
        f.write(f"{Nobs}\n")
        f.write("\n")
        f.write("\n")
        while t < t1:
            f.write(f"{t} {Ndata}\n")
            f.write(f"{x_r} {y_r} {z_r}\n")
            f.write(f"{x_o} {y_o} {z_o}\n")
            # Dummy flux
            f.write(f"10 10000 100\n")
            f.write(f"\n")
            t += dt_d
    # Observation file ========================================================
