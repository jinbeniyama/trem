#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Common functions to prepare TPM and to handle/plot/utilize the results.

Read document carefully!
JPL/HORIZONS: https://ssd.jpl.nasa.gov/horizons/manual.html
astroquery  : https://astroquery.readthedocs.io/en/latest/jplhorizons/jplhorizons.html

- About location
location is like '(MPCcode)' (e.g., 381), 
                 '(MPCcode)@(399, which means Earth)' (e.g., 381@391)
The latter part is to specify major body 
    (Earth '391', Sun '10', Moon '301', solar system barycenter '0' or 'ssb' etc.)
The former part is to specify point on the major body (body center 500 etc.)

J.B. notes that we have to set '@10', or 'None' (not '@0' or 'ssb')
since  the center of Ecliptic coordinate system, 
which is used in Marco's TPM code, is the Sun (not solar system barycenter).


location refers to the coordinate center for the ephemeris, which has 
slightly different physical interpretations depending on the query type: 
  - observer (ephemerides) queries
  - observer location vectors queries
  - (coordinate origin for vectors elements queries)
  - (relative body for orbital elements)

The default value is None for all queries, which corresponds to
  - Earth body center for observer (ephemerides) queries 
    (i.e., location=500, 500@391)
  - Sun body center for orbital elements and vectors queries 
    (i.e., location=@10, 500@10)


2024-05-15 J.B. confirmed that the generated ephem file with @10 
           matched that in DAFEED better than that with @0. 
           J.B. confirmed that the generated obs file with light time correction
           matched that in DAFEED with an accuracy of 0.01 s.
"""
import numpy as np
from astroquery.jplhorizons import Horizons


mycolor = [
    "#AD002D", "#1e50a2", "#69821b", "#f055f0", "#afafb0", 
    "#0095b9", "#89c3eb", "#ec6800", "cyan", "gold", "magenta"
    ] 
mycolor = mycolor*500
mymark = ["o", "^", "s", "D", "*", "v", "<", ">", "h", "H"]
mymark = mymark*500



def make_ephemfile(asteroid, df, out, warmuptime_day=30):
    """
    Make ephem file for TPM.
    Light time correction is unnecessary.

    Parameters
    ----------
    asteroid : str
        name of asteroid
    df : pandas.DataFrame
        jd, wavelength, flux, fluxerr, code, and cflag
    out : str
        output ephem filename
    warmuptime_day : float, optional
        time for warming up in day
    """
    # Light-time correction is needless! (private com. with Marco DELBO, July 19 2024)

    # Sort by jd
    jds = sorted(set(df["jd"]))

    # Make chunks. A chunk is composed of obs. with shorter arc
    # This is just to speed up the process.
    # th_arc (50 days) 
    #   i.e., if observation gap is larger than 70 days, 
    #         they are regarded as different chunks
    # Note: th_arc should be larger than margin0 + margin1, otherwise
    #       ephmeris file could be not always ascending order
    th_arc = 50
    jds_chunk = []
    
    # First day
    d0, d1 = jds[0], jds[0]
    for idx_d, d in enumerate(jds):
        if d - d0 < th_arc:
            d1 = d
            # For final day
            if idx_d == len(jds)-1:
                # Save first and last days as a chunk
                ch = (d0, d1)
                jds_chunk.append(ch)
        else:
            # Save first and last days as a chunk
            ch = (d0, d1)
            jds_chunk.append(ch)

            # Update first day
            d0, d1 = d, d

            # For final day
            if idx_d == len(jds)-1:
                # Save first and last days as an independent chunk
                ch = (d, d)
                jds_chunk.append(ch)


    # Make ephem file widh jds_chunk ==========================================
    # Use 1 day time step. Each obs. is linked by interpolation.
    # Just to speed up,  all jds are not used.
    # t0 of chunk 1: (minimum jd of chunk 1) - 50, -50 is margin to reach thermal equilibrium.
    #                Note that this may be enough for slow rotators.
    # t1 of chunk 1: (maximum jd of chunk 1) + 1, 1 is just in case
    #                And add 1 flag to move to next chunk (i.e., t0 of chunk 2)
    # t0 of chunk 2: (minimum jd of chunk 2) - 50, -50 is margin to reach thermal equilibrium.
    # t1 of chunk 2: (maximum jd of chunk 2) + 1, 1 is just in case
    #                And add 1 flag to move to next chunk (i.e., t0 of chunk 3)
    # ......
    # 
    # t0 of final chunk : (minimum jd of final chunk) - 50, -50 is margin to reach thermal equilibrium.
    # t1 of final chunk : (maximum jd of final chunk) + 1, 1 is necessary to run a tpm I don't know why.

    # Note: margin0 of 30 corresopnds to 60 rotations for object wiht rotation period of 12 hr.
    margin0, margin1 = warmuptime_day, 1
    with open(out, "w") as f_eph:
        for idx_ch, ch in enumerate(jds_chunk):
            d0, d1 = ch
            # Add margins
            t0 = d0 - margin0
            t1 = d1 + margin1

            d_list = np.arange(t0, t1, 1)
            for idx, d in enumerate(d_list):
                # Location of @10 means Sun body center (=None) for vectors queries
                # (not solar system barycenter, @0, @ssb).
                S = Horizons(location="@10", id=asteroid, epochs=d)

                # 2024-10-21 for test
                #S = Horizons(location="@0", id=asteroid, epochs=d)

                vec = S.vectors(refplane="ecliptic")
                # Vector from the Sun to asteroid (Sun -> ast)
                x_S, y_S, z_S = vec["x"][0], vec["y"][0], vec["z"][0]
                # Vector from asteroid to the Sun (ast -> Sun)
                x_S, y_S, z_S = -x_S, -y_S, -z_S
                # Last raw should end with 1 to skip (or end) the calculations.
                if (idx == len(d_list)-1):
                    f_eph.write(f"{d} {x_S} {y_S} {z_S} 1\n")
                else:
                    f_eph.write(f"{d} {x_S} {y_S} {z_S}\n")
     

def make_obsfile(asteroid, df, out, lccor=False, rmnegativeflux=False):
    """
    Parameters
    ----------
    asteroid : str
        name of asteroid
    df : pandas.DataFrame
        jd, wavelength, flux, fluxerr, code, and cflag
    out : str
        output obs filename
    lccor : bool, optional
        wheather perform light-time correction
    """
    # Number of data block or epoch
    N_epoch = len(set(df["jd"]))

    # Sort by jd
    jds = sorted(set(df["jd"]))
    # 0,   1,  2,  3,       4?,     5?,      6?,       7?,       8,        9
    # W1, W2, W3, W4, IRAS 12, IRAS 25, IRAS 60, IRAS 100, Akari 9, Akari 18
    cflag_registered = ["0", "1", "2", "3", "8", "9", ]

    # Make obs file ===========================================================
    with open(out, "w") as f_obs:
        f_obs.write(f"{N_epoch}\n")
        f_obs.write(f"\n")
        for jd in jds:
            # Extract 1 data point. (i.e., merge data points if it is spectrum)
            df_samejd = df[df["jd"]==jd]
            df_samejd = df_samejd.reset_index(drop=True)
            # Number of observations
            N_data = len(df_samejd)

            # Common parameters for all data (i.e., spectrum)
            code, cflag   = df_samejd["code"][0], df_samejd["cflag"][0]
            if str(cflag) in cflag_registered:
                pass
            else:
                cflag = ""

            
            # Location of @10 means Sun body center (=None) for vectors queries
            # (not solar system barycenter, @0, @ssb).
            S = Horizons(location="@10", id=asteroid, epochs=jd)
            # 2024-10-21 for test
            #S = Horizons(location="@0", id=asteroid, epochs=jd)

            vec = S.vectors(refplane="ecliptic")
            # Vector from the Sun to asteroid (Sun -> ast)
            x_S, y_S, z_S = vec["x"][0], vec["y"][0], vec["z"][0]
            # Vector from asteroid to the Sun (ast -> Sun)
            x_S, y_S, z_S = -x_S, -y_S, -z_S

            # Location of code means observer center
            E = Horizons(location=code, id=asteroid, epochs=jd)
            vec = E.vectors(refplane="ecliptic")
            eph = E.ephemerides()

            # Vector from the Earth to asteroid (Earth -> ast)
            x_E, y_E, z_E = vec["x"][0], vec["y"][0], vec["z"][0]
            # Vector from asteroid to the Earth (ast -> Earth)
            x_E, y_E, z_E = -x_E, -y_E, -z_E
            
            if lccor:
                # Lighttime correction
                # Geocentric distance in au
                delta = eph["delta"][0]
                # c: speed of light in m/s
                # au: astronomical unit in m
                c_au_s = c.value/au.value
                c_au_day = c_au_s*24.*3600.
                ltcor = delta/c_au_day
                jd_ltcor = jd - ltcor
            else:
                print("    Do not perform light time correction!")
                jd_ltcor = jd
            
            # Count number of negative flux
            if rmnegativeflux:
                N_nega = 0
                for idx, row in df_samejd.iterrows():
                    flux = row["flux"]
                    if flux < 0:
                        N_nega += 1
                N_data -= N_nega

            # Write 
            #   jd, Ndata
            #   Sun coordinate
            #   Earth coordinate
            f_obs.write(f"{jd_ltcor} {N_data}\n")
            f_obs.write(f"{x_S} {y_S} {z_S}\n")
            f_obs.write(f"{x_E} {y_E} {z_E}\n")
            
            for idx, row in df_samejd.iterrows():
                w             = row["wavelength"]
                flux, fluxerr = row["flux"], row["fluxerr"]
                # Do not use negative flux
                if (rmnegativeflux) and (flux < 0):
                    continue
                f_obs.write(f" {w} {flux} {fluxerr} {cflag}\n")

            f_obs.write("\n")
