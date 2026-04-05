#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Common functions for dual-component TPM.
"""
from trem.emittance.common_emittance import calc_chi2, calc_chi2_numpy
import numpy as np

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
    
    # 2. This is faster =======================================================
    df_blend = df1.copy()
    df_blend["f_model"] = (
        alpha*df1["scalefactor"]**2*df1["f_model"] + 
        (1-alpha)*df2["scalefactor"]**2*df2["f_model"])
    # This is a dummy
    df_blend["scalefactor"] = 1
    # 2. This is faster. =======================================================

    return df_blend


def blend_flux_numpy(f1, s1, f2, s2, alpha):
    """
    Blend fluxes if df1 and df2 with alpha.
      F_blended = alpha*F1*s1^2 + (1 - alpha)*F2*s2^2,
    where s1 and s2 are scale factors.
    """
    return alpha * (s1**2) * f1 + (1 - alpha) * (s2**2) * f2


#def search_regolith_abundance(df1, df2, alpha_list, chi2_min=10000, minonly=False):
#    """
#    Search regolith abundance alpha which minimize chi2.
#
#    Parameters
#    ----------
#    df1 : pandas.DataFrame
#        dataframe with TI of regolith (i.e., low TI)
#    df2 : pandas.DataFrame
#        dataframe with TI of regolith (i.e., low TI)
#    alpha_list : array-like
#        list of regolith abundance
#    chi2_min : float
#        initial chi2 minimum
#    minonly : bool
#        return minimum chi2 and corresponding alpha (i.e., fit by alpha)
#    sf : bool
#        introduce scale parameters per epoch (only for spectra)
#
#    Returns
#    -------
#    alpha_arr : float
#        array of regolith abundance 
#    chi2_arr : float
#        array of chi2
#    """
#    alpha_arr, chi2_arr = [], []
#    
#    f1 = df1["f_model"].to_numpy()
#    s1 = df1["scalefactor"].to_numpy()
#    f2 = df2["f_model"].to_numpy()
#    s2 = df2["scalefactor"].to_numpy()
#    f_obs = df1["f_obs"].to_numpy()
#    ferr_obs = df1["ferr_obs"].to_numpy()
#
#    for a in alpha_list:
#        # Blend flux as 
#        #   F = alpha*F_regolith*s1^2 + (1-alpha)*F_rock*s2^2,
#        # where s1 and s2 are scale factors.
#        ## This is slow
#        #df_blend = blend_flux(df1, df2, a)
#
#        ## This is faster
#        f_blend = blend_flux_numpy(f1, s1, f2, s2, a)
#
#        # Calculate chi2 of blended flux
#        ## This is slow
#        #chi2 = calc_chi2(df_blend)
#        ## This is faster
#        ## Set global scale factor to 1 (scale factors are alraeady introduced!)
#        chi2 = calc_chi2_numpy(f_obs, f_blend, ferr_obs, 1)
#
#        if minonly:
#            if chi2 < chi2_min:
#                alpha_arr = [a]
#                chi2_arr = [chi2]
#                chi2_min = chi2
#            else:
#                pass
#        else:
#            alpha_arr.append(a)
#            chi2_arr.append(chi2)
#
#    return alpha_arr, chi2_arr

# Fast
def search_regolith_abundance(df1, df2, alpha_list, chi2_min=10000, minonly=False):
    """
    Search regolith abundance alpha which minimize chi2.
    Vectorized implementation (same algorithm, faster).
    """

    # convert to numpy arrays
    f1 = df1["f_model"].to_numpy()
    s1 = df1["scalefactor"].to_numpy()
    f2 = df2["f_model"].to_numpy()
    s2 = df2["scalefactor"].to_numpy()
    f_obs = df1["f_obs"].to_numpy()
    ferr_obs = df1["ferr_obs"].to_numpy()

    alpha_arr = []
    chi2_arr = []

    # Precompute scaled fluxes
    rego = f1 * s1**2
    rock = f2 * s2**2

    alpha = np.asarray(alpha_list)

    # Vectorized blend flux: F = rock + alpha * (rego - rock)
    delta = rego - rock
    f_blend = rock + alpha[:, None] * delta  # shape: (N_alpha, N_data)

    # Vectorized chi2
    diff = (f_obs - f_blend)**2 / ferr_obs**2
    chi2_all = np.sum(diff, axis=1)

    if minonly:
        idx = np.argmin(chi2_all)
        if chi2_all[idx] < chi2_min:
            alpha_arr = [alpha[idx]]
            chi2_arr = [chi2_all[idx]]
    else:
        alpha_arr = alpha.tolist()
        chi2_arr = chi2_all.tolist()

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
