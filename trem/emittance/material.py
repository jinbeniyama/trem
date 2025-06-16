#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions related to material properties.
Many functions are copy from the public code. 
- regolith_conductivity.dvrc by Andrew J. Ryan
Notations are modified by J.B.

Notations
---------
k_s (wo/suffi  x)   : solid thermal conductivity
k_rad               : radiative thermal conductivity
Phi                 : microporosity (inside particles)
phi                 : macroporosity (between particles)
"""
from astropy import constants as const
SB_const = const.sigma_sb.value


# Functions to derive k_s (solid thermal conductivity) from T (temperature) ==
def k_s_glass(T):
    """
    Glass solid conductivity based on a modeled curve (Ratcliff, 1963 model, S-Type Jaygo).
    
    Parameters
    ----------
    T : float
        Temperature in Kelvin.
        
    Returns
    -------
    k_s_glass : float
        Thermal conductivity in W/m·K.
    """
    # Fit to modeled glass curve based on Ratcliff, 1963 model. S-Type Jaygo. 
    k_s_glass = 0.3319 + (4.216e-3)*T - (5.1e-6)*T**2
    return k_s_glass


def k_s_glass5(T):
    """
    Glass solid conductivity for 5mm diameter glass based on Ratcliff, 1963 model (M-Type Jaygo).
    
    Parameters
    ----------
    T : float
        Temperature in Kelvin.
        
    Returns
    -------
    k_s_glass5 : float
        Thermal conductivity in W/m·K.
    """
    # Fit to modeled glass curve based on Ratcliff, 1963 model. M-Type Jaygo (5 mm diameter)
    k_s_glass5 = 0.2852 + (3.175e-3)*T - (2.7e-6)*T**2
    return k_s_glass5


def k_s_glass10(T):
    """
    Glass solid conductivity for 10mm diameter glass based on Ratcliff, 1963 model (M-Type Jaygo).
    
    Parameters
    ----------
    T : float
        Temperature in Kelvin.
        
    Returns
    -------
    k_s_glass10 : float
        Thermal conductivity in W/m·K.
    """
    # Fit to modeled glass curve based on Ratcliff, 1963 model. M-Type Jaygo (1 cm diameter)
    k_s_glass10 = 0.2296 + (4.7e-3)*T - (6.05e-6)*T**2
    return k_s_glass10


def k_s_basalt(T):
    """
    Basalt solid conductivity from Desai et al., 1974 (used by Sakatani et al., 2016, 2018).
    
    Parameters
    ----------
    T : float
        Temperature in Kelvin.
        
    Returns
    -------
    k_s_basalt : float
        Thermal conductivity in W/m·K.
    """
    k_s_basalt = (-9.53e-4)*T + 2.4
    return k_s_basalt


def k_s_tagish(T):
    """
    Tagish lake sim conductivity estimated from Cold Bokkeveld (CM2), Opeil et al 2010.
    
    Parameters
    ----------
    T : float
        Temperature in Kelvin.
        
    Returns
    -------
    k_s_tagish : float
        Thermal conductivity in W/m·K.
    """
    # simple fit, valid from 100-300 K, can may be extrapolated to higher temps
    k_s_tagish = 0.26 + 0.0013*T
    return k_s_tagish


def k_s_CK(T):
    """
    NWA 5515 (CK4) conductivity from Opeil et al., 2010.
    
    Parameters
    ----------
    T : float
        Temperature in Kelvin.
        
    Returns
    -------
    k_s_CK : float
        Thermal conductivity in W/m·K.
    """
    # NWA 5515 (CK4)
    k_s_CK = 1.26 + 0.0011*T
    return k_s_CK


def k_s_ottawa(T):
    """
    Ottawa sand conductivity based on a model from Ottawa sand data.
    
    Parameters
    ----------
    T : float
        Temperature in Kelvin.
        
    Returns
    -------
    float
        Thermal conductivity in W/m·K.
    """
    k_s_ottawa = (3880.8)*(T**-1.0667)
    return k_s_ottawa


def k_s_glass_sakatani(T):
    """
    Glass solid conductivity from Sakatani et al., 2017.
    
    Parameters
    ----------
    T : float
        Temperature in Kelvin.
        
    Returns
    -------
    k_s_glass_sakatani : float
        Thermal conductivity in W/m·K.
    """
    k_s_glass_sakatani = (8.50e-4)*T + 0.855
    return k_s_glass_sakatani


def k_s_graphite(T):
    """
    Graphite conductivity based on van Antwerpen et al. (2012).
    Input temperature is in Celsius, so it is converted to Kelvin inside the function.
    
    Parameters
    ----------
    T : float
        Temperature in Kelvin.
        
    Returns
    -------
    k_s_graphite : float
        Thermal conductivity in W/m·K.
    """
    # From van Antwerpen et al 2012 figure 24 and 25
    T_celsius = T - 273.15  # Convert to Celsius
    k_s_graphite = 186.021 - (39.5408e-2)*T_celsius + (4.8852e-4)*T_celsius**2 - (2.91e-7)*T_celsius**3 + (6.6e-11)*T_celsius**4
    return k_s_graphite
# Functions to derive k_s_solid (thermal conductivity) from T (temperature) ==


# Dictionary of material properties with their corresponding conductivity functions
material_properties = {
    "glass": 
        {"poissons": 0.22, "youngs": 63.e9, "rho_bulk": 2495.0, "k_s_func": k_s_glass},
    "glass5": 
        {"poissons": 0.22, "youngs": 65.e9, "rho_bulk": 2495.0, "k_s_func": k_s_glass5},
    "glass10": 
        {"poissons": 0.22, "youngs": 63.e9, "rho_bulk": 2570.0, "k_s_func": k_s_glass10},
    "basalt": 
        {"poissons": 0.25, "youngs": 78.e9, "rho_bulk": 3040.0, "k_s_func": k_s_basalt},
    "sakatani_basalt": 
        {"poissons": 0.25, "youngs": 78.e9, "rho_bulk": 2900.0, "k_s_func": lambda T: (-9.53e-4)*T + 2.40},
    "tagish": 
        {"poissons": 0.269, "youngs": 5.625e9, "rho_bulk": 2270.0, "k_s_func": k_s_tagish},
    "CK": 
        {"poissons": 0.269, "youngs": 5.625e9, "rho_bulk": 2675.0, "k_s_func": k_s_CK},
    "sakatani_glass": 
        {"poissons": 0.22, "youngs": 55.1e9, "rho_bulk": 2480.0, "k_s_func": k_s_glass_sakatani},
    # Same with above
    "FGB": 
        {"poissons": 0.22, "youngs": 55.1e9, "rho_bulk": 2480.0, "k_s_func": k_s_glass_sakatani},
    "EMB": 
        {"poissons": 0.22, "youngs": 55.1e9, "rho_bulk": 2600.0, "k_s_func": lambda T: (5.10e-4)*T + 1.406},
    # Valid at 200-300 K, based on porosity fit from Opeil 2012, assuming LL meteorite mean 8.2% porosity
    "itokawa": 
        {"poissons": 0.25, "youngs": 78.e9, "rho_bulk": 3220.0, "k_s_func": lambda T: 1.40},
    # Valid at 200-300 K, based on porosity fit from Opeil 2012, assuming LL meteorite mean 8.2% porosity
    "eros": 
        {"poissons": 0.25, "youngs": 78.e9, "rho_bulk": 3220.0, "k_s_func": lambda T: 1.40},
    "ottawa": 
        {"poissons": 0.08, "youngs": 1.10e11, "rho_bulk": 2650.0, "k_s_func": k_s_ottawa}
}



def get_material_properties(sample, T):
    """
    Get material properties based on the sample name and temperature.
    
    Parameters
    ----------
    sample : str
        Name of the material sample (e.g., 'glass', 'basalt', etc.)
    T : float
        Temperature (K)
    
    Returns
    -------
    dict
        Dictionary containing material properties:
        'poissons' : Poisson's ratio,
        'youngs' : Young's Modulus (Pa),
        'rho_bulk' : Bulk density (kg/m³),
        'k_s' : Thermal conductivity (W/m·K).
    """
    sample = sample.lower()
    if sample in material_properties:
        material = material_properties[sample]
        k_s = material["k_s_func"](T)
        material["k_s"] = k_s
        return material
    else:
        raise ValueError(f"Sample '{sample}' is not supported!")


# Radiative component of thermal conductivity =================================
def calc_k_rad_GB(e, diam, T, phi, e1):
    """
    Calculate the radiative component of regolith conductivity based on the Gundlach & Blum model (2013).

    Parameters
    ----------
    e : float
        Emissivity of the material
    diam : float
        Diameter of the particle (m?)
    T : float
        Temperature (K)
    phi : float
        macroporosity of the material
    e1 : float
        scaling factor for radiative heat transfer?.

    Returns
    -------
    k_rad_GB : float
        Radiative conductivity (W/m·K).
    """
    k_rad_GB = 8.0 * SB_const * e * (T**3) * e1 * (phi / (1.0 - phi)) * (diam / 2.0)
    return k_rad_GB


def calc_k_rad_sakatani(e, diam, T, phi, zeta):
    """
    Calculate the radiative thermal conductivity based on the Sakatani et al. (2017) model.

    Parameters
    ----------
    e : float
        Emissivity of the material
    diam : float
        Diameter of particle (m)
    T : float
        Temperature (K)
    phi : float
        macroporosity of the material
    zeta : float
        Scaling factor related to particle size?

    Returns
    -------
    k_rad_sakatani : float
        Radiative thermal conductivity (W/m·K).
    """

    # Values for zeta from Sakatani et al., 2017
    # 710-1000 um: 0.7-1.0
    # 355-500 um: 1.1-1.9
    # 180-250 um: 1.2-1.7
    # 90-106 um: 1.8-2.6
    # 53-63 um: 2.5-4.0
    # emb powder 5 um: ~15
    k_rad_sakatani = 8.0 * (e / (2.0 - e)) * SB_const * zeta * ((phi / (1.0 - phi))**(1.0 / 3.0)) * (diam / 2.0) * T**3
    return k_rad_sakatani
# Radiative component of thermal conductivity =================================


def c_p_ordinary_chondrite(type_meteorite, T_typical):
    """
    Heat capacity of ordinary chondrites from Mache et al. (2019), MPS, 54, 2729.
    Valid from 70 to 300 K.

    Parameters
    ----------
    type_meteorite : str
        H, L, LL
    T_typical : float
        temperature in K

    Return
    ------
    c_p : float
        heat capacity in J/kg/K
    """
    
    assert 70 < T_typical < 300, "Invalid temperature."
     
    if type_meteorite == "H":
        A, B, C, D = (-213.4, 5.659, -11.131e-3, 8.697e-6)
    elif type_meteorite == "L":
        A, B, C, D = (-206.4, 5.507, -9.684e-3, 6.258e-6)
    elif type_meteorite == "LL":
        A, B, C, D = (-200.9, 5.394, -8.815e-3, 4.787e-6)
    else:
        assert False, "Not implemeted other than H, L, LL."
    # Heat capacity
    c_p = A + B*T_typical + C*T_typical**2 + D*T_typical**3
    return c_p
