#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Common functions for dual-component TPM.
"""
from trem.emittance.common_emittance import calc_chi2


def blend_flux(df1, df2, alpha):
    """
    Blend fluxes if df1 and df2 with alpha.
      F_blended = alpha*F1*s1^2 + (1 - alpha)*F2*s2^2,
    where s1 and s2 are scale factors.

    Parameters
    ----------
    df1 : pandas.DataFrame
        dataframe with TI of regolith (i.e., low TI)
    df2 : pandas.DataFrame
        dataframe with TI of regolith (i.e., low TI)
    alpha : float
        regolith abundance

    Return
    ------
    df_blend : pandas.DataFrame
        dataframe with blended results
    """
    # Sanity check
    assert len(df1) == len(df2), "Check if input dfs are the same dimension!"
    
    # 1. This is safe but slow, because of the for loop. ======================
    # Number of data points
    #N_data = len(df1)
    #epoch_list, jd_list, w_list = [], [], []
    #f_obs_list, ferr_obs_list, f_model_list = [], [], []
    #
    #for n in range(N_data):
    #    epoch1    = df1.loc[n, "epoch"]
    #    jd1       = df1.loc[n, "jd"]
    #    w1        = df1.loc[n, "w"]
    #    f_obs1    = df1.loc[n, "f_obs"]
    #    ferr_obs1 = df1.loc[n, "ferr_obs"]
    #    f_model1  = df1.loc[n, "f_model"]
    #    s1        = df1.loc[n, "scalefactor"]

    #    epoch2    = df2.loc[n, "epoch"]
    #    jd2       = df2.loc[n, "jd"]
    #    w2        = df2.loc[n, "w"]
    #    f_obs2    = df2.loc[n, "f_obs"]
    #    ferr_obs2 = df2.loc[n, "ferr_obs"]
    #    f_model2  = df2.loc[n, "f_model"]
    #    s2        = df2.loc[n, "scalefactor"]

    #    # Sanity checks
    #    assert jd1 == jd2, "Check if input dfs are made from the TPMs with same obs file."
    #    assert w1 == w2, "Check if input dfs are made from the TPMs with same obs file."
    #    assert f_obs1 == f_obs2, "Check if input dfs are made from the TPMs with same obs file."
    #    assert ferr_obs1 == ferr_obs2, "Check if input dfs are made from the TPMs with same obs file."

    #    # Calculate blended flux
    #    f_blend = alpha*s1**2*f_model1 + (1 - alpha)*s2**2*f_model2

    #    # Save info.
    #    epoch_list.append(epoch1)
    #    jd_list.append(epoch1)
    #    w_list.append(w1)
    #    f_obs_list.append(f_obs1)
    #    ferr_obs_list.append(ferr_obs1)
    #    f_model_list.append(f_blend)
    #    
    ## DataFrame
    #df_blend = pd.DataFrame({
    #    "epoch": epoch_list, 
    #    "jd": jd_list, 
    #    "w": w_list,
    #    "f_obs": f_obs_list, 
    #    "ferr_obs": ferr_obs_list, 
    #    "f_model": f_model_list, 
    #    })
    ## This is a dummy
    #df_blend["scalefactor"] = 1
    # 1. This is safe but slow, because of the for loop. ======================

    # 2. This is faster =======================================================
    df_blend = df1.copy()
    df_blend["f_model"] = (
        alpha*df1["scalefactor"]**2*df1["f_model"] + 
        (1-alpha)*df2["scalefactor"]**2*df2["f_model"])
    # This is a dummy
    df_blend["scalefactor"] = 1
    # 2. This is faster. =======================================================

    return df_blend


def search_regolith_abundance(df1, df2, alpha_list, chi2_min=10000, minonly=False):
    """
    Search regolith abundance alpha which minimize chi2.

    Parameters
    ----------
    df1 : pandas.DataFrame
        dataframe with TI of regolith (i.e., low TI)
    df2 : pandas.DataFrame
        dataframe with TI of regolith (i.e., low TI)
    alpha_list : array-like
        list of regolith abundance
    chi2_min : float
        initial chi2 minimum
    minonly : bool
        return minimum chi2 and corresponding alpha (i.e., fit by alpha)
    sf : bool
        introduce scale parameters per epoch (only for spectra)

    Returns
    -------
    alpha_arr : float
        array of regolith abundance 
    chi2_arr : float
        array of chi2
    """
    alpha_arr, chi2_arr = [], []
    
    for a in alpha_list:
        # Blend flux as 
        #   F = alpha*F_regolith*s1^2 + (1-alpha)*F_rock*s2^2,
        # where s1 and s2 are scale factors.
        df_blend = blend_flux(df1, df2, a)
        # Calculate chi2 of blended flux
        chi2 = calc_chi2(df_blend)

        if minonly:
            if chi2 < chi2_min:
                alpha_arr = [a]
                chi2_arr = [chi2]
                chi2_min = chi2
            else:
                pass
        else:
            alpha_arr.append(a)
            chi2_arr.append(chi2)

    return alpha_arr, chi2_arr


def calc_C_coord(phi):
    """
    Calculate coordination number (see Sakatani+2018)

    Parameter
    ---------
    phi : float
        macroscopic porosity
        
    Return
    ------
    C_coord : float
        coordination bumber C 
    """
    f = 0.07318 + 2.193 * phi
    C_coord = 2.812 * (1 - phi)**(-1./3.) / (f**2 * (1 + f**2))
    return C_coord


def kappa_Sakatani2018(kappa, D_p, phi, r_c, xi):
    """
    Calculate bulk thermal conductivity 
    with model in Sakatani+2018, Icarus, 309, 13.

    Parameters
    ----------
    kappa : float
        thermal conductivity of solid material
    D_p : float
        particle diameter
    phi : float
        macroscopic porosity
    r_c : float
        radius of the contact area between the spheres
    xi : float
        degree of reduction of the thermal conductance at the contacts 
        owring to the microscopic roughness of the particle surfaces

    Return
    ------
    kappa_bulk : float
        bulk thermal conductivity
    """
    # Particle radius
    R_p = D_p/2.

    # Ratio of 
    # [the effective distance of radiative heat transfer in the voids between particles] 
    # to [the geometric size of the voids]
    # used in Cambioni+2021
    ## 0.68 + 7.6e-5 / D_p
    ## ???
    ## How to calculate r_c?? Assume Young's modulus, Poisson's ratio?
    r_c = xxxx

    # Calculate Coordination number C with phi
    C_coord = calc_C_coord(phi)

    # Equation (8) in Sakatani+2018
    kappa_bulk = 4 / np.pi**2 * kappa * (1 - phi) * C_coord * xi * r_c/R_p
    return kappa_bulk


def calc_Phi(TI_rock, c_p, rho_s):
    """
    Calculate microscopic porosity Phi.

    Parameters
    ----------
    TI_rock : float
        thermal inertia of rock 
    c_p : float
        heat capacity
    rho_s : float
        grain density

    Return
    ------
    Phi : float
        microscopic porosity
    """
    # TI_rock**2/(c_p rho_s (1-Phi)) = 0.11(1-Phi)/Phi
    # Define C = (TI_rock**2)/(0.11 c_p rho_s) + 2 and solve the equation
    # (See note on Cambioni+2021 by JB)
    C = TI_rock**2/(0.11*c_p*rho_s) + 2
    # Solution satisfying 0 < Phi < 1
    Phi = (C - (C**2 - 4)**0.5) / 2
    return Phi
    

def calc_TI_th(TI_rock):
    """
    Calculate threshold of thermal inertia of regolith and rock
    (gamma_c in Cambioni+2021).
    gamma_c is defined as thermal inertia when D_p = l_s,
    where D_p is particle diameter and l_s is thermal skin depth.

    Parameters
    ----------
    TI_rock : float
        thermal inertia of rock

    Return
    ------
    TI_th : float
        threshold of thermal inertia of regolith and rock
    """

    # Fixed parameters in Cambioni+2021 =======================================
    # Macroporosity of Bennu 
    phi = 0.4
    # Grain density of CM meteorites in kg/m^3
    rho_s = 2920
    # Heat capacity for meteorite CM2 Cold Bokkeveld at OTES spot's 
    # mean diurnal temperature (from Figure 3 in Opeil+2020, MPS)
    # TODO: what is mean diurnal temperature? Read by eye?
    c_p = 999
    # Degree of reduction of the thermal conductance at the contacts 
    # owring to the microscopic roughness of the particle surfaces
    # (from Sakatani+2018, Icarus ?)
    xi = 0.12
    # Fixed parameters in Cambioni+2021 =======================================


    # 1. Derive microscipic porosity (microporosity) Phi
    #    with TI_rock, c_p (heat capacity), and rho_s (grain density)
    Phi = calc_Phi(TI_rock, c_p, rho_s)

    # 2. Derive thermal conductivity kappa
    kappa = TI_rock**2/(c_p * rho_s * (1 - Phi))

    # Density considering macroporosity phi and microscopic porosity Phi
    rho = rho_s * (1 - Phi) * (1 - phi)

    # 3.Derive regolith bulk thermal conductivity, kappa_p,
    # with Standard Regolith Model in Sakatani+2018
    # TODO: This kappa is dependent on D_p?
    kappa_p_list = kappa_Sakatani2018(kappa)
    #       Find k_p where D_p = l_s?
    kappa_p = 999


    # 4. Finally derive threshold of thermal inertia
    TI_th = (kappa_p * c_p * rho)**0.5
    return TI_th
