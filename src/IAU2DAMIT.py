#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Conert IAU spin to DAMIT spin.

Originally this code was written by Joel Beccarelli. 
Original matlab code was provided by Josef Durech.
The code is validated with a MATLAB code written by Saverio Cambioni.

Please see the documentation on DAMIT.
https://damit.cuni.cz/projects/damit/pages/documentation
For the sake of clarity, we use t1 rather than t on the DAMIT.

Note on the time:
    DAMIT: UTC
    IAU  : Barycentric Dynamical Time (Temps Dynamique Barycentrique, TDB)

Note on the consistency:
J.B. guesses that 
the inconsistency between values in DAMIT or Table 2 in Mottola+2020
and the outputs of this code are from the fact that some digits 
have been rounded in the references. 
For Ceres, for instance, if we use 
  (alpha, delta = 290.4, 63.48)
rather than 
  (alpha, delta = 290, 63), # values in DAMIT IAUspin.txt
we can get consistent phi0.
For Leucus in Mottola+2020, if we use
  (alpha, delta = 247.6, 58.15)
rather than 
  (alpha, delta = 248, 58), # values in their Table 2
we can get consistent phi0.


Examples
--------
Validation with J.B.'s favorite asteroid (433) Eros

>>> python3 IAU2DAMIT.py --obj Eros

Calculate DAMIT spin parameters with arbitrary inputs

>>> python IAU2DAMIT.py --alpha 290 --delta 63 --W0 247.3 --W1 952.1529 --t0 2451545.00 --t1 2434407.00

# Phobos from Archinal+2011
>>> python IAU2DAMIT.py --alpha 317.68 --delta 52.90 --W0 35.06 --W1 1128.8445850 --t0 2451545.00 --t1 2451545.00
"""
from argparse import ArgumentParser as ap
import numpy as np
from numpy.linalg import inv
from astropy.time import Time


# Useful functions ============================================================
def tdb2utc(jd_tdb: float) -> float:
    """Convert Julian Date in TDB (Barycentric Dynamical Time) to UTC.

    Parameters
    ----------
    jd_tdb : float
        Julian Date in TDB

    Returns
    -------
    t_utc: float
        Julian Date in UTC
    """
    t_tdb = Time(jd_tdb, format="jd", scale="tdb")
    t_utc = t_tdb.utc
    t_utc = t_utc.jd
    return t_utc


def fetch_param(obj: str) -> dict:
    """Return parameters for validation.

    Parameter
    ---------
    obj : str or int
        object name or number

    Return
    ------
    params : dictionary
        dict with parameters
    """
    # Capital just to judge
    obj = obj.upper() 
    # Ref: DAMIT
    if obj == "EROS":
        alpha = 11
        delta = 17
        W0 = 326.1   
        W1 = 1639.389365 
        t0 = 2451545.00
        t1 = 2451545  
        phi0  = 32.64
    # Ref: DAMIT
    # This does not reproduce the values on DAMIT.
    # (delta phi0 = 2.8106 deg)
    elif obj == "CERES":
        alpha = 290
        delta = 63
        
        # Value to get consistent phi0
        #alpha = 290.4
        #delta = 63.48

        W0    = 247.3   
        W1    = 952.152878
        t0    = 2451545.00
        t1    = 2434407
        phi0  = 0
    # Ref: Used in Cambioni+2021, Nature
    elif obj == "BENNU":
        # I didn't check these values
        alpha = 86.6388
        delta = -65.1086
        W0 = 89.64
        # What is W1 from SPICE Kernel?
        P_hr = 4.296003
        W1 = 360 / P_hr * 24
        t0 = 2451545.00
        t1 = 2451545.00
        # Value in spinDelbo2019.05.27.txt
        phi0  = 122.31
        assert False, "Check the code"
    # Ref: DAMIT
    # This does not reproduce the values on DAMIT.
    # (delta phi0 = 0.3527 deg)
    elif obj == "VESTA":
        alpha = 314
        delta = 42
        W0 = 94.1
        W1 = 1617.334280
        t0 = 2451545.00
        t1 = 2433574.0
        phi0  = 0
    # Ref: DAMIT
    # This does not reproduce the values on DAMIT.
    # (delta phi0 = 0.1433 deg)
    elif obj == "220622":
        alpha = 54
        delta = -13
        W0 = 18.0
        W1 = 1397.625961
        t0 = 2451545.00
        t1 = 2456908
        phi0  = 0
    # Ref: Carry+2010
    # This does not reproduce the values on DAMIT.
    # (delta phi0 ~ 2.3 deg)
    elif obj == "LUTETIA":
        alpha = 52
        delta = 12
        W0    = 94
        # = 360 / 8.168270 * 24
        W1    = 360. / 8.168270 * 24.
        t0    = 2451545.00
        t1    = 2444822.35116
        phi0  = 0
    # Ref: Mottola+2020, PSJ, 1, 73.
    # (delta phi0 ~ 0.7 deg)
    elif obj == "LEUCUS":
        alpha = 248
        delta = 58

        # Value to get consistent phi0
        alpha = 247.6
        delta = 58.15

        W0    = 60.014   
        W1    = 360/445.683*24
        t0    = 2451545.00
        t1    = 2456378
        phi0  = -76.129%360

    else:
        raise NotImplementedError(f"Parameters for {obj} are not prepared yet.")

    params = {
        # See DAMIT doc for details
        # https://damit.cuni.cz/projects/damit/pages/documentation
        "target": obj,
        "alpha": alpha,
        "delta": delta,
        "W0": W0,     
        "W1": W1,
        "t0": t0,   
        "t1": t1,
        "phi0": phi0
        }
    return params


def Rx(r): 
    """Equatoreal frame, x axis, right
    """
    return np.array([
    [1, 0, 0],
    [0, np.cos(r), np.sin(r)],
    [0, -np.sin(r), np.cos(r)]
])

def Ry(r):
    """y axis, right
    """
    return np.array([
        [np.cos(r), 0, -np.sin(r)],
        [0, 1, 0],
        [np.sin(r), 0, np.cos(r)]
])

def Rz(r):
    """z axis, right
    """
    return np.array([
        [np.cos(r), np.sin(r), 0],
        [-np.sin(r),  np.cos(r), 0],
        [0, 0, 1]
])


def cart2sph(x, y, z):
    """Convert Cartesian coordinates (x, y, z) to spherical coordinates.
    
    Parameters
    ----------
    x : float 
        X coordinate
    y : float 
        Y coordinate
    z : float 
        Z coordinate.

    Returns
    -------
    azi, elev, r
    """
    r = np.sqrt(x**2 + y**2 + z**2)
    # angle in xy-plane from x-axis
    azi = np.arctan2(y, x)             
    # angle from xy-plane toward z-axis
    elev = np.arcsin(z / r)          
    return azi, elev, r
# Useful functions ============================================================


# Main ========================================================================
def IAU2DAMIT(
    alpha: float, delta: float, W0: float, W1: float, t0: float, t1: float):
    """Main function to convert IAU spin to DAMIT one. 

    Parameters
    ----------
    alpha : float
        equatorial longitude of the pole in deg
    delta : float
        equatorial latitude of the pole in deg
    W0 : float
        the position of the prime meridian (positive x) at the time of t0 in deg 
    W1 : float
        rotation rate in deg/day (IAU definition)
    t0 : float
        initial epoch for which the position of the prime meridian is W0
    t1 : float
        epoch for which the position of the prime meridian is W1
        (t on DAMIT)

    Return
    ------
    params : dictionary
        lam (ecliptic longitude), beta (ecliptic latitude), 
        P_hr (sidereal rotation period in hour), and phi0 (rotation phase at t0)
    """
    # Axial tilt for the Earth
    e0_deg = 23 + 26 / 60 + 21.4119 / 3600
    e0_rad = np.radians(e0_deg)
    R_eps = Rx(e0_rad)
    
    pole_body = [0, 0, 1]
    x_body = [1, 0, 0]
    
    # Rotation period in hour
    P_hr = 1. / (W1 / 360.) * 24.
    
    pole_eq = [np.cos(np.radians(delta))*np.cos(np.radians(alpha)), 
               np.cos(np.radians(delta))*np.sin(np.radians(alpha)),
               np.sin(np.radians(delta))]
    pole_ecl = R_eps@pole_eq
    
    lam, beta, r = cart2sph(pole_ecl[0], pole_ecl[1], pole_ecl[2])
    
    # lam and beta are already radian
    if lam < 0:
        lam = lam + (2 * np.pi)

    lam_deg = np.degrees(lam)
    beta_deg = np.degrees(beta)

    # This doesn't change that much
    #W0 = (W0 + W1*(t1 - t0))%360
    #t0 = tdb2utc(t0)

    # W: position of the prime meridian at the time t in deg
    # Added by J.B.
    W = (W0 + W1*(t1 - t0))%360
    
    M_ = inv(
        Rz(np.radians(W)) @ Rx(np.radians(90-delta)) @ Rz(np.radians(90+alpha)))
    x_eq = M_@x_body

    R_eps_nega    = Rx(-e0_rad)
    R_lambda = Rz(-lam)
    R_beta    = Ry(-(np.pi/2.-beta))

    M__ = R_eps_nega @ R_lambda @ R_beta
    x_ecl = inv(M__)@x_eq
    az, el, radius = cart2sph(x_ecl[0], x_ecl[1], x_ecl[2])
    phi0 = np.degrees(az)

    params = {
            "lam": lam_deg,
            "beta": beta_deg,
            "P_hr": P_hr,
            "phi0": phi0%360.,
        }
    return params
# Main ========================================================================


if __name__ == "__main__":
    parser = ap(description="Convert IAU spin to DAMIT spin")
    parser.add_argument(
        "--obj", type=str, default=None,
        help="Object name for validation (Eros)")
    parser.add_argument(
        "--alpha", type=float, default=290,
        help="Equatorial long of north pole (in IAU)")
    parser.add_argument(
        "--delta", type=float, default=63,
        help="Equatorial lat of north pole (in IAU)")
    parser.add_argument(
        "--W1", type=float, default=952.152878,
        help="rotation rate in deg/day (in IAU)")
    parser.add_argument(
        "--t0", type=float, default=2451545.00,
        help="t0 (time in IAU)")
    parser.add_argument(
        "--W0", type=float, default=247.3,
        help="phase W0 of the PM in J2000 (in IAU)")
    parser.add_argument(
        "--t1", type=float, default=2434407.00,
        help="t1 (time in DAMIT spin file)")
    args = parser.parse_args()
    

    # Validation
    if args.obj:
        params = fetch_param(args.obj)
        alpha, delta = params["alpha"], params["delta"]
        W0, W1 = params["W0"], params["W1"]
        t0, t1 = params["t0"], params["t1"]
        phi0 = params["phi0"]

        params_DAMIT = IAU2DAMIT(
            alpha, delta,
            W0, W1, 
            t0, t1)

        print()
        print(f"Calculate DAMIT spin parameter")
        print(f"Validation for {args.obj}")
        print()

        print("   Inputs")
        print(f"    (alpha, delta) = ({alpha}, {delta}) deg")
        print(f"    W0             = {W0:.4f} deg")
        print(f"    W1             = {W1:.4f} deg/day")
        print(f"    t0             = {t0:.4f}")
        print(f"    t1             = {t1:.4f}")
        print(f"    Correct phi0   = {phi0:.4f} deg")
        print()

        print("   Calculated")
        print(f"    rotP_hr     = {params_DAMIT['P_hr']:.4f} hr")
        print(f"    phi0        = {params_DAMIT['phi0']:.4f} deg")
        print(f"    (lam, beta) = ({params_DAMIT['lam']:.4f}, {params_DAMIT['beta']:.4f}) deg")
        print()
   
    else:
        params_DAMIT = IAU2DAMIT(
            args.alpha, args.delta,
            args.W0, args.W1,
            args.t0, args.t1)

        print()
        print(f"Calculate DAMIT spin parameter")
        print()
       
        print("   Inputs")
        print(f"    (alpha, delta) = ({args.alpha}, {args.delta}) deg")
        print(f"    W0             = {args.W0:.4f} deg")
        print(f"    W1             = {args.W1:.4f} deg/day")
        print(f"    t0             = {args.t0:.4f}")
        print(f"    t1             = {args.t1:.4f}")
        print()

        print("   Calculated")
        print(f"    rotP_hr   = {params_DAMIT['P_hr']:.4f} hr")
        print(f"    phi0      = {params_DAMIT['phi0']:.4f} deg")
        print(f"    lam       = {params_DAMIT['lam']:.4f} deg")
        print(f"    beta      = {params_DAMIT['beta']:.4f} deg")
        print()
        pass
