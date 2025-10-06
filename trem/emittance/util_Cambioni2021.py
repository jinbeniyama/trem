#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions mainly from Cambioni et al., 2021, Nature.
Many functions are copy from the public code. 
- regolith_conductivity.dvrc by Andrew J. Ryan
Some functions (e.g., TI_range, particle_size) are not used in the determination of 
the thermal inertia cutoff.
Notations are modified by J.B.

Notations
---------
k_s (wo/suffi  x)   : solid thermal conductivity
k_rad               : radiative thermal conductivity
Phi                 : microporosity (inside particles)
phi                 : macroporosity (between particles)
rho_s               : grain density
rho_b               : bulk density, i.e., (1.0 - Phi) * rho_s
rho_e               : effective bulk density, i.e., (1 - phi)(1 - Phi) * rho_s
"""
from argparse import ArgumentParser as ap
import numpy as np
from astropy import constants as const
import matplotlib.pyplot as plt
from trem.emittance.material import (
    get_material_properties, calc_k_rad_GB, calc_k_rad_sakatani, 
    c_p_ordinary_chondrite)

# Constants
# Stefan-Boltzmann constant (W/m^2/K^4)
SB_const = const.sigma_sb.value


def calc_prop(desired_TI, specific_heat=750.0, model="avg", rho_s=2920.0, rotP_hr=4.296057, phi=0.50):
    """
    Computes the relevant thermal conductivity, density, and related properties
    based on the desired thermal inertia (TI).

    Parameters
    ----------
    desired_TI : float
        Target thermal inertia value.
    specific_heat : float, optional
        Specific heat value (default: 750.0).
    model : str, optional
        Model for converting between porosity and conductivity. Options are
        "avg", "henke", and "flynn" (default: "avg").
    rho_s : float, optional
        Grain density, default is 2920.0 (average CM from Macke et al., 2011).
    phi : float, optional
        macroscopic porosity

    Returns
    -------
    dict
        A dictionary containing:
        - k : Thermal conductivity.
        - rho : Density.
        - phi : Porosity.
        - sigma_t : Surface roughness temperature fluctuation scaling factor.
        - skin_depth : Thermal skin depth.
    """
    # Rotational period in seconds
    rotP_s = rotP_hr * 60.0 * 60.0  

    # microporosity values
    Phi = np.linspace(0.002, 0.98, 1000)

    # Thermal conductivity models
    flynn_k = 0.11 * (1.0 - Phi) / Phi
    henke_k = 4.3 * np.exp(-Phi / 0.08)
    avg_k = (flynn_k + henke_k) / 2.0

    print(f"  flynn_k, henke_k, avg_k = {flynn_k[0]}, {henke_k[0]}, {avg_k[0]}")

    if model == "henke":
        avg_k = henke_k
    elif model == "flynn":
        avg_k = flynn_k

    # Bulk density
    rho_b = rho_s * (1 - Phi)
    # Effective bulk density
    rho_e = rho_b * (1 - phi)

    # Thermal inertia calculation
    # Note: rho_b is used in the original code (Cambioni+2021)
    TI = np.sqrt(rho_b * avg_k * specific_heat)

    # Find porosity corresponding to the desired TI
    if desired_TI < TI.min() or desired_TI > TI.max():
        raise ValueError("Desired TI is out of the calculated range.")

    out_Phi = np.interp(desired_TI, TI[::-1], Phi[::-1])
    out_k = np.interp(out_Phi, Phi, avg_k)
    out_rho_b = rho_s * (1.0 - out_Phi)

    # Surface roughness temperature fluctuation scaling factor
    # Serpentine value from Grott et al., 2019
    ks = 2.95  
    sigma_t = (np.pi / 4000.0) * (out_k / ks) * (10.0e6)

    # Thermal skin depth
    skin_depth = np.sqrt(out_k * rotP_s / (out_rho_b * specific_heat * np.pi))


    return {
        "k": out_k,
        "rho_b": out_rho_b,
        "Phi": out_Phi,
        "sigma_t": sigma_t,
        "skin_depth": skin_depth,
    }


def Nc_Suzuki(phi):
    """
    Calculate the coordination number (Nc).

    Parameters
    ----------
    phi : float
        macroporosity of the material (between 0 and 1).

    Returns
    -------
    Nc : float
        The coordination number Nc.
    """
    # Used by Sakatani 2017; 2018 for calculating coordination number. 
    f = 0.07318 + 2.193*phi - 3.357*phi**2 + 3.194*phi**3
    # coordination number, Suzuki et al 1981
    Nc = (2.812*(1.0 - phi)**(-1./3.))/((f**2.)*(1.0 + f**2.)) 
    return Nc


def make_Fmech2(depth, rho_e, Rs, phi, planet):
    """
    Calculate the mechanical force due to depth, density, radius, and porosity for a given planet.
    
    Parameters
    ----------
    depth : float
        Depth of the material in m.
    rho_e : float
        Effective bulk density of the material in kg/m^3
        rho_e = (1-phi)(1-Phi)rho_s, where rho_s is material density.
    Rs : float
        Radius of the spherical object in m.
    phi : float
        Macroscopic porosity
    planet : str
        Name of the planet (e.g., 'Earth', 'Mars', 'Moon', etc.).
    
    Returns
    -------
    float
        Mechanical force (Newtons).
    """
    # Set gravitational acceleration based on the planet
    planet_gravity = {
        "earth": 9.80665,
        "mars": 3.711,
        "moon": 1.622,
        "itokawa": 0.086e-3,
        "bennu": 75.0e-6,
        "eros": 0.003,
        "ryugu": 1.6e-4
    }

    # Convert planet name to lowercase for case-insensitive comparison
    planet = planet.lower()
    
    # Get gravitational acceleration, default to Earth if planet not supported
    g = planet_gravity.get(planet, 9.80665)
    
    # Calculate the mechanical force
    Fmech = 2.0 * np.pi * rho_e * g * depth * (Rs**2) / (np.sqrt(6.0) * (1.0 - phi))
    
    return Fmech


def rc_jkr(Rloc, Fmech, Gamma, v, E):
    """
    Calculate the relative contact radius using the JKR (Johnson-Kendall-Roberts) model.
    
    Parameters
    ----------
    Rloc : float
        Radius of the particle or spherical object (meters).
    Fmech : float
        Force due to the overlying mechanical load (Newtons).
    Gamma : float
        Surface energy (J/m²).
    v : float
        Poisson's ratio (dimensionless).
    E : float
        Young's Modulus (Pascals).
    
    Returns
    -------
    float
        Relative contact radius (dimensionless).
    """
    Rstar = Rloc / 2.0
    F_con_jkr = (Fmech + 3 * np.pi * Gamma * Rstar + np.sqrt(6 * np.pi * Gamma * Rstar * Fmech + (3 * np.pi * Gamma * Rstar)**2))
    R_con_jkr = (3 * Rstar * (1 - v**2) * F_con_jkr / (2 * E))**(1 / 3.0)
    return R_con_jkr


def keff(sphere_diam, depth, distance, phi, T, P, emiss, rho, gas_type, 
         fk=1.0, rc_factor=1.0, sample="glass", xi=0.60, X=0.41, zeta=1.0, 
         e1=1.34, k_const=None, rc_fixed=None, rc_fixed_ratio=None, planet="Earth", 
         verbosity=0, Nc_const=None, contact_angle_const=None, new_fk=0, zetaxi=None, surfenergy=0.032):
    """
    Calculate the effective thermal conductivity (keff) of the material.

    Parameters
    ----------
    sphere_diam : float
        Diameter of the sphere in meters. 
        This is the characteristic size of the particles or grains.
    depth : float
        Depth in meters. 
        Used for calculating the mechanical force (lithostatic pressure).
    distance : float
        Distance in meters. 
        Not used in the calculation but could represent the distance between particles.
    phi : float
        macroporosity of the material (dimensionless). 
        The fraction of void spaces within the material (between 0 and 1).
    T : float
        Temperature in Kelvin. 
        Temperature at which thermal conductivity is evaluated.
    P : float
        Pressure in Pascals. 
        This parameter is provided for completeness but not used in the current calculation.
    emiss : float
        Emissivity of the material (dimensionless). 
        Describes how efficiently the material radiates energy.
    rho : float
        Density of the material in kg/m³. 
        The mass per unit volume of the material.
    gas_type : str
        Type of gas (optional). 
        This parameter is not used in the current calculation but may be included for future extensions.
    fk : float, optional
        Scaling factor for radiative heat transfer (default is 1.0).
        A correction factor for radiative heat transfer based on material properties.
    rc_factor : float, optional
        Scaling factor for the contact radius (default is 1.0). 
        A factor to adjust the calculated contact radius.
    sample : str, optional
        Type of material ('glass' or 'basalt', default is 'glass').
        Determines the properties of the material used for calculations.
    xi : float, optional
        Shape factor for particles (default is 0.60). 
        A factor to account for the shape of particles in the material.
    X : float, optional
        Gas-specific parameter (default is 0.41). 
        A parameter related to the specific gas properties (not used here).
    zeta : float, optional
        Scaling factor for radiative conductivity (default is 1.0).
        A factor for adjusting radiative conductivity based on porosity.
    e1 : float, optional
        Scaling factor for radiative heat transfer (default is 1.34).
        A parameter used in the Gundlach & Blum radiative heat transfer model.
    k_const : float, optional
        Fixed thermal conductivity value (default is None).
        If provided, this value will override the calculated thermal conductivity.
    rc_fixed : float, optional
        Fixed contact radius (default is None).
        If provided, this value will override the calculated contact radius.
    rc_fixed_ratio : float, optional
        Contact radius ratio (default is None).
        If provided, this will be used to adjust the contact radius.
    planet : str, optional
        Name of the planet (default is 'Earth'). 
        Used in calculating the mechanical force (lithostatic pressure).
    verbosity : int, optional
        Level of verbosity for detailed output (default is 0).
        If set to 1, it will print detailed information about the calculation.
    Nc_const : float, optional
        Fixed coordination number (default is None).
        If provided, this value will override the calculated coordination number.
    contact_angle_const : float, optional
        Fixed contact angle (default is None).
        A fixed contact angle, if known, to be used in calculations.
    new_fk : int, optional
        Flag for using a new fk calculation method (default is 0).
        If set to 1, it will use a new method for the fk parameter calculation.
    zetaxi : float, optional
        Relationship between zeta and xi (default is None).
        Used if the zeta and xi parameters are related and need to be adjusted.
    surfenergy : float, optional
        0.032 is for glass

    Returns
    -------
    out : array-like
    out["k_s_sakatani"] = k_s_sakatani
    out["k_rad_sakatani"] = k_rad_sakatani
    out["k_s_GB"] = k_s_GB
    out["k_rad_GB_rad"] = k_rad_GB
    out["rc"] = rc
    out["Nc"] = Nc
    out["rho"] = rho
    out["fk_predicted"] = fk_predicted
    out["zeta"] = zeta
    out["xi"] = xi
    """
 
    # Coordination number with macroporosity, phi
    Nc = Nc_Suzuki(phi)

    # Override coordination number if provided
    if Nc_const is not None:
        Nc = Nc_const  

    # Get material properties
    props = get_material_properties(sample, T)
    poissons = props["poissons"]
    youngs = props["youngs"]
    k_s_func = props["k_s_func"]
    rho_b = props["rho_bulk"]

    k_s = k_s_func(T)
    
    # Consider macroporosity
    rho_e = rho_b * (1.0 - phi)

    # Calculate the mechanical force (lithostatic pressure)
    Fmech2 = make_Fmech2(depth, rho_e, sphere_diam / 2.0, phi, planet)
    rc_sakatani = rc_jkr(sphere_diam / 2.0, Fmech2, surfenergy, poissons, youngs) * rc_factor
    
    if rc_fixed:
        rc_sakatani = rc_fixed
    if rc_fixed_ratio:
        rc_sakatani = (sphere_diam/2.0) * rc_fixed_ratio

    gamma = rc_sakatani / sphere_diam
    
    # The same as Cambioni+2021
    if zetaxi:
        xi = 0.12
    

    # Solid conductivity calculations
    Hs = ((4.0 * np.pi / 3.0) ** (1.0 / 3.0)) * (sphere_diam / 2.0) * k_s
    Hctotal = Nc * (2.0 / np.pi) * k_s * rc_sakatani * xi
    H = 1.0 / ((1.0 / Hs) + (1.0 / Hctotal))
    k_s_sakatani = (2.0 / np.pi) * (1.0 - phi) * H / (sphere_diam / 2.0)

    
    # Gundlach & Blum Solid conductivity
    surfenergySiO = (6.67e-5) * T  # Equation for surface energy
    f1 = 5.18e-2
    f2 = 5.26
    buffer = f1 * np.exp(f2 * (1.0 - phi))
    k_s_GB = (((9.0 * np.pi / 4.0) * ((1.0 - poissons ** 2.0) / youngs) * (surfenergySiO / (sphere_diam / 2.0))) ** (1.0 / 3.0)) * buffer * X * k_s

    # Radiation: Non-isothermal correction factor (fk)
    Lambda_s = k_s / (4.0 * sphere_diam * SB_const * T ** 3)
    if verbosity == 1:
        print(f"Isothermal correction Lambda: {1.0 / Lambda_s}")

    if new_fk == 0:
        # Parameters from Van Antwerpen et al., 2012
        a1 = 0.0841 * emiss ** 2 - 0.307 * emiss - 0.1737
        a2 = 0.6094 * emiss + 0.1401
        a3 = 0.5738 * emiss ** (-0.2755)
        a4 = 0.0835 * emiss ** 2 - 0.0368 * emiss + 1.0017
    else:
        # Values from Ryan et al., 2020
        a1 = -0.568
        a2 = 0.912
        a3 = 0.765
        a4 = 1.035

    fk_predicted = a1 * np.arctan(a2 * (1 / Lambda_s) ** a3) + a4

    # TODO: Ensure fk <= 1?
    # This seems always 1
    fk_predicted = np.where(fk_predicted > 1.0, 1.0, fk_predicted)
    #print(fk_predicted)

    if verbosity == 1:
        print(f"fk predicted: {fk_predicted}")

    if zeta is None:
        zeta = 1.0
    if zetaxi == 1:
        zeta = 0.68 + (7.6e-5) / sphere_diam  # From Wada et al., 2018


    # Radiative conductivity from Sakatani model
    k_rad_sakatani = calc_k_rad_sakatani(emiss, sphere_diam, T, phi, zeta)  
    
    # Radiative conductivity from Gundlach & Blum model
    k_rad_GB = calc_k_rad_GB(emiss, sphere_diam, T, phi, e1)  

    rc = rc_sakatani

    out = {}
    out["k_s_sakatani"] = k_s_sakatani
    out["k_rad_sakatani"] = k_rad_sakatani
    out["k_s_GB"] = k_s_GB
    out["k_rad_GB"] = k_rad_GB
    out["rc"] = rc
    out["Nc"] = Nc
    out["rho_e"] = rho_e
    out["fk_predicted"] = fk_predicted
    out["zeta"] = zeta
    out["xi"] = xi
    return out


def calc_TIth(TI_rock, T_typical, obj, phi):
    """Calculate threshold of thermal inertia of regolith and rock, gamma_c.

    gamma_c is defined as thermal inertia when D_p = l_s,
    where D_p is particle diameter and l_s is thermal skin depth.

    Parameters
    ----------
    TI_rock : float
        thermal inertia of rock in tiu
    T_typical : float
        typical temperature in K
    obj : str
        target name (e.g., Bennu, Eros)
    phi: float
        Macroscopic porosity

    Returns
    -------
    Dth : float
        particle diameter
    TIth : float
        threshold of thermal inertia of regolith and rock
    """

    # Fixed parameters in Cambioni+2021 =======================================
    # Non-isothermal correction factor
    fk = 1
    # Relationship between and k (conductivity) and phi 
    # "flynn" is used in Cambioni+2021
    model = "flynn"
    # Fixed parameters in Cambioni+2021 =======================================


    if obj == "Bennu":
        # Parameters in Cambioni+2021
        # Emissivity (fixed to 0.95 for Bennu in Cambioni+2021)
        emiss = 0.95

        # Rotation period in hr
        rotP_hr = 4.296057
        rotP_s = rotP_hr*3600.

        # Macroporosity
        # 0.15, 0.40 (nominal), and 0.60 are used.

        # Grain density of CM meteorites in kg/m^3
        rho_s = 2920

        # Heat capacity for meteorite CM2 Cold Bokkeveld at OTES spot's 
        # mean diurnal temperature (from Figure 3 in Opeil+2020, MPS)
        # OTES spot's diurnal temperature is from ~220 up to 320 
        # (Supp. material of Rotizits+2020, Sci. Adv., 6, eabc3699., abc3699_sm.pdf)
        # as a function of (T_typical)
        # So constant value cannot reproduce the 
        # same result as Cambioni+2021
        # T ~ 300 K 
        c_p = 750.0


    elif obj =="Eros":
        # Emissivity (The same as TPM in Beniyama+)
        emiss = 0.90

        # Rotation period in hr
        rotP_hr = 5.27
        rotP_s = rotP_hr*3600.

        # Macroporosity
        # Wilkison+2002, Icarus, 155, 94
        # best estimated to be 20%
        #phi = 0.2

        # Grain density of ordinary meteorites in kg/m^3
        # From table 1 of Macke et al. 2019, MPS, 54, 2729.
        rho_s = 3600

        # At Eros's mean "GLOBAL?" diurnal temperature (from Figure 4 of Macke+2019, MPS, 54, 2729.)
        meteorite_type = "H"
        c_p = c_p_ordinary_chondrite(meteorite_type, T_typical)
    else:
        assert False, "Not implemented."

    # Degree of reduction of the thermal conductance at the contacts 
    # owring to the microscopic roughness of the particle surfaces
    # (from Sakatani+2018, Icarus ?)
    xi = 0.12

    # Particle diameter array from 100 microns to 15 cm
    # in m
    #D_arr = np.linspace(0.100e-3, 0.150e-2, 150)
    D_arr = np.linspace(100e-6, 150e-3, 10000)

    # Calculate k_m (conductivity) and rho (material density of rock fragments) based on TIrock
    result = calc_prop(
        TI_rock, model=model, specific_heat=c_p, rho_s=rho_s, rotP_hr=rotP_hr, phi=phi)  
    # rho_b: Bulk density
    k_m, rho_b = result["k"], result["rho_b"]
    # rho_e: Effective bulk density
    rho_e = rho_b * (1-phi)

    # Create lookup table of regolith conductivity vs particle size
    out_keff = keff(
        D_arr, 0.01, 10.0, phi, T_typical, 1.0e-10, emiss, rho_e, 
        "N2", sample="tagish", planet=obj, k_const=k_m, new_fk=1, zetaxi=1, surfenergy=0.032) 
    
    k_s_sakatani   = out_keff["k_s_sakatani"]
    k_rad_sakatani = out_keff["k_rad_sakatani"]
    fk_predicted = out_keff["fk_predicted"]

     # Calculate corrected conductivity based on fk
    if fk > 0:
        k_out = k_s_sakatani + k_rad_sakatani * fk_predicted
    else:
        k_out = k_s_sakatani + k_rad_sakatani

    # Calculate skin depth using regolith conductivity vs diameter
    skin = np.sqrt(k_out * rotP_s / (rho_e * c_p * np.pi))

    # Find the intersection point of skin depth and particle size curves
    mindifpos = np.argmin(np.abs(D_arr - skin))  
    # Find the closest value of D to skin depth
    Dth = D_arr[mindifpos]
    k_out = k_out[mindifpos]
    TIth = np.sqrt(k_out * c_p * rho_e)

    return Dth, TIth



if __name__ == "__main__":
    parser = ap(description="Test to check the functions.")
    parser.add_argument(
        "action", type=str, default="TIth",
        help="Plot TI rock vs. TI th")
    parser.add_argument(
        "--obj", type=str, default="Bennu",
        help="Target object")
    parser.add_argument(
        "--T_typical", type=float, default=300.,
        help="Typical temperature")
    args = parser.parse_args()

    
    # Plot TIrock vs. TIth
    if args.action == "TIth":
        obj = args.obj
        T_typical = args.T_typical
        TIrock_list = np.arange(25, 2500, 25)

        print(f"  Plot TIrock vs. TIth")
        print(f"  Object: {obj}, T={T_typical} [K]")
        print(f"  Temperature: {T_typical}")

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_axes([0.15, 0.15, 0.7, 0.7])
        ax.set_xlabel("TI of rock [tiu]")
        ax.set_ylabel("TI cutoff [tiu]")
        ax.set_xscale("log")
        ax.set_title(f"{obj}, T={T_typical} [K]")
        ax.set_ylim([0, 200])

        phi_list = [0.15, 0.40, 0.60]
        ls_list = ["solid", "dashed", "dotted"]
        col_list = ["black", "red", "blue"]
        for idx, phi in enumerate(phi_list):
            TIth_list = []
            for TI_rock in TIrock_list:
                _, TIth = calc_TIth(TI_rock, T_typical, obj, phi)
                TIth_list.append(TIth)

            ax.plot(
                TIrock_list, TIth_list, ls=ls_list[idx], color=col_list[idx], label=f"$\phi={phi}$")

        ax.legend(loc="lower right")
        plt.savefig(f"TIrock_vs_TIth_{obj}.jpg")
    
    # For future update
    elif args.action == "aaa":
        pass
    else:
        print("  Not implemented : {args.action}")
