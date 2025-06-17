#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Common functions to prepare TPM and to handle/plot/utilize the results (emittance).
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from trem.common import mycolor, mymark


def extract_flux(f0, fixscale=False):
    """
    Extract thermal flux from output of TPM.

    Parameters
    ----------
    f : str
        output file of TPM
    fixscale : bool, optional
        whether fix the scale factor (i.e., trust shape model)

    Return
    ------
    df : pandas.DataFrame
        dataframe with extracted fluxes etc.
    """

    epoch_list, jd_list, w_list = [], [], []
    f_obs_list, ferr_obs_list, f_model_list = [], [], []
    with open (f0, "r") as f:
        f = f.readlines()
        for l in f:
            # l is as below:
            # f> 052     2454199.2663881136   000 18.000        1.698843   0.113268    1.441020    0.000000      1.1789
            #    epoch   JD                   n   wavelength    f_obs      ferr_obs    f_model     ferr_model    f_obs/f_model
            if l[0:2] == "f>":
                # Extract info
                l = l.split()
                epoch, jd, n       = l[1], l[2], l[3]
                w, f_obs, ferr_obs = l[4], l[5], l[6]
                f_model            = l[7]
                epoch_list.append(float(epoch))
                jd_list.append(float(jd))
                w_list.append(float(w))
                f_obs_list.append(float(f_obs))
                ferr_obs_list.append(float(ferr_obs))
                f_model_list.append(float(f_model))
            elif l[0:2] == "r>":
                # Extract scale factor
                # r>     100.0   0.0   0.0  0.00  0.12 1.15125549     19.364     21.576  1.000   2259.587    4.18547
                l = l.split()
                # Fix scale factor to 1
                if fixscale:
                    scalefactor = 1
                # Use scale factor in TPM output
                else:
                    scalefactor = float(l[6])

                # This is useful to check intputs of the TPM calculation
                # Haple angle in deg (t.Bar in the output) 
                TI = float(l[1])
                Htheta = float(l[2])
                A = float(l[5])

            else:
                continue

        # DataFrame
        df = pd.DataFrame({
            "epoch": epoch_list,
            "jd": jd_list,
            "w": w_list,
            "f_obs": f_obs_list,
            "ferr_obs": ferr_obs_list,
            "f_model": f_model_list,
            })
        df["scalefactor"] = scalefactor
        df["TI"] = TI
        df["Htheta"] = Htheta
        df["A"] = A
        #print(f"Scafe factor = {scalefactor}")
    return df


def plot_flux(dfs, fixscale=False, y1range=None, out=None):
    """
    Plot obs and model fluxes.
    Note: model fluxes of input files should be the same.

    Parameters
    ----------
    dfs : array-like
        preprocessed dataframes with fluxes
    out : str, optional
        output filename
    """
    fig = plt.figure(figsize=(24, 16))
    # idx vs. flux (obs. and model, normalized by obs.)
    ax1 = fig.add_axes([0.05, 0.70, 0.90, 0.25])
    # residual
    ax2 = fig.add_axes([0.05, 0.45, 0.90, 0.25])
    # chi2 component 
    ax3 = fig.add_axes([0.05, 0.20, 0.90, 0.25])

    ax3.set_xlabel("Index")
    ax1.set_ylabel("Normalized flux")
    if fixscale:
        ax2.set_ylabel("(f_obs - f_model)/ferr_obs")
    else:
        ax2.set_ylabel("(f_obs - sf**2*f_model)/ferr_obs")
    ax3.set_ylabel("Chi2 component")

    for idx_df, df in enumerate(dfs):
        col = mycolor[idx_df]
        idx_epoch = 0
        chi2 = 0

        for idx_row, row in df.iterrows():

            epoch    = row["epoch"]
            f_obs    = float(row["f_obs"])
            ferr_obs = float(row["ferr_obs"])
            f_model  = float(row["f_model"])
            s        = float(row["scalefactor"])

            # Fit with the size
            if fixscale:
                pass
            else:
                f_model = f_model*s**2

            # Set marker
            if idx_row == 0:
                pass
            else:
                # Change marker for different epoch
                if epoch != epoch0:
                    idx_epoch += 1

            col = mycolor[idx_epoch]
            mark = mymark[idx_epoch]

            # f_obs and s**2 f_model
            # ax.set_ylabel("Flux [Jy]")
            #ax.scatter(idx_row, f_obs, color=col, marker=mark)
            #ax.errorbar(idx_row, f_obs, ferr_obs, color=col)
            #ax.scatter(idx_row, s**2*f_model, color=col, marker=mark)
            
            # Scale factor is already taken into account in f_model above
            # f_obs - s**2 f_model/ferr_obs
            if f_obs > 0:
                res = (f_obs - f_model)/ferr_obs
                chi2 += res**2
            print(f"  {idx_row+1}/{len(df)} f_obs, f_model, res = {f_obs:.5f}+-{ferr_obs:.5f}, {f_model:.5f}, {res:.5f}")

            if idx_row == (len(df)-1):
                label = f"Idx {idx_df}: chi2 = {chi2:.3f}"
                label_obs = f"Observations (normalized by themselves, i.e., 1)"
                label_model = f"Models (normalized by observations)"
            else:
                label = None
                label_obs = None
                label_model = None

            # Observation
            ax1.errorbar(idx_row, f_obs/f_obs, ferr_obs/f_obs, color=col, ms=0, marker=mark, label=label_obs)
            # Model
            ax1.scatter(idx_row, f_model/f_obs, color="black", s=5, marker="o", label=label_model)

            # Residual
            ax2.scatter(idx_row, res, color=col, marker=mark, label=label)
            # Chi2 component
            ax3.scatter(idx_row, res**2, color=col, marker=mark, label=label)

            # Update previous epoch
            epoch0 = epoch
    
    if y1range:
        ax1.set_ylim([y1range[0], y1range[1]])
    ymin, ymax = ax2.get_ylim()
    if abs(ymin) > abs(ymax):
        ax2.set_ylim([ymin, -ymin])
    else:
        ax2.set_ylim([-ymax, ymax])

    ax1.legend(ncol=3)
    ax2.legend(ncol=3)
    if out:
        plt.savefig(out)


def crater2Htheta(c_angle, c_density):
    """
    Calculate Haple roughness parameter with 
    semiaperture angle of craters and crater surface density.
    See Hapke 1984, Icarus.
    # Note: This function is no longer useful since
    #       I realised that Htheta was written in the output of TPM!

    Parameters
    ----------
    c_angle : float
        semiaperture angle of craters in degree
    c_density : float
        crater surface density (0 to 1)

    Return
    ------
    Htheta : float
        Haple roughness parameter in degree
        
    """
    # TODO: Calculate theta_bar by myself
    # How to calculate theta_bar?
    Htheta = c_density**0.5*np.tan(np.radians(c_angle))

    # This is temporary one.
    # From Hung+2022, PSJ, Table 2
    # Convert c_angle to theta_bar
    # theta_bar is 12.0 for (gamma, rhoC) = (50.0, 0.50), not 22.0
    # (i.e., Hung+2022 is correct.)
    d = {0.0:0.0, 30.0:3.9, 40.0:12.6, 41.0:16.5, 50.0:12.0, 
         60.0:26.7, 70.0:27.3, 88.0:35.8, 89.0:46.8, 90.0:55.4}
    Htheta = d[c_angle]

    return Htheta
 

def extract_bestparam(df, key_chi2, params):
    """
    Search and return best fit parameters.

    Parameters
    ----------
    df : pandas.DataFrame
        input dataframe
    key_chi2 : str
        key of chi-squared
    params : array-like
        key of parameters of interest 

    Returns
    -------
    chi2_min : float
        minimum chi-squared
    best_params : array-like
        best fit parameters and chi2 minimum
    """
    # Extract minimum chi2 and its index
    idx_min = df[key_chi2].idxmin()
    chi2_min = df.loc[idx_min, key_chi2]

    best_params = []
    for key in params:
        best_param = df.loc[idx_min, key]
        best_params.append(best_param)
    return chi2_min, best_params


def calc_chi2(df):
    """
    Calculate chi2 with F_obs, Ferr_obs, F_model, and scale factor.
    Note: 
    The output chi2 could be different from that in the output of TPM.
    This is because runtpm ignores negative fluxes. J.B. thinks negative 
    fluxes are still informative if they are with large errorbars.

    Parameter
    ---------
    df : pandas.DataFrame
        dataframe with fluxes

    Return
    ------
    chi2 : float
        calculated chi-square
    """
    df["diff"] = (df["f_obs"] - df["scalefactor"]**2*df["f_model"])**2/df["ferr_obs"]**2
    chi2 = np.sum(df["diff"])
    return chi2


def calc_confidence_chi2(paper, chi2_min, dof, n, reduce):
    """
    Calculate condifence levels.
    There are several ways to estimate confidence levels,
    and errors of physical properties using the confidence levels.

    1. P14
       Polishook 2014 Icarus, 241, 79, Cambioni+2021, Nature, 
       chi2 = chi2_min + n * (2/nu)**0.5 (for n-sigma)

    2. V17
       Vokrouhlicky+2017, AJ, Hanus+2018, A&A, Durech+2018 A&A, etc.
       chi2 = chi2_min * (1 + n*(2/nu)**0.5)
            = chi2_min + chi2_min*n*(2/nu)**0.5 (for n-sigma)
       Note: 
           Cambioni+2019 as well with n=1, Equatioin (4) is a typo,
           private com. with Saverio on 2025-01-09
    
    Both of these two expressions are approximations.
    (chi2_min is not the unity in reality)
    Since normally chi2_min > 1, 1. gives smaller confidence levels and smaller 
    uncertainties of physical proproperties compared to 2.
    When chi2_min ~ 1, the results are almost the same. 
    See confidence_level.ipynb for comparison.

    Parameters
    ----------
    paper : str
        P14 or V17
    chi2_min : float
        minimum chi-squared
    dof : int
        degrees of freedom
    n : int
        return n-sigma interval
    reduce : bool
        whether reduced chi-squared or not

    Return
    ------
    chi2_sigma : 
        chi-squared boundary with n-sigma confidence level
    """

    if reduce:
        chi2_sigma = n*np.sqrt(2*dof)/dof
    else:
        chi2_sigma = n*np.sqrt(2*dof)

    # 1. Polishook 2014, Cambioni+2021, 
    if paper == "P14":
        pass
    # 2. Vokrouhlicky+2017, AJ, Hanus+2018, A&A, Durech+2018 A&A, etc.
    elif paper == "V17":
        chi2_sigma = chi2_sigma*chi2_min

    return chi2_sigma


def introduce_var_scalefactor(df, key_t="jd", sf0=0.80, sf1=1.2, sfstep=0.01):
    """
    Introduce variable scale factors per observation.

    Parameters
    ----------
    df : pandas.DataFrame
        input dataframe
    key_t : str
        keyword for observation time
    sf0 : float
        minimum scale factor        
    sf1 : float
        maximum scale factor        
    sfstep : float
        step of scale factor

    Return
    ------
    df1 : pandas.DataFrame
        output dataframe with new scale factors
    """

    # Updated dataframe
    df1_list = []

    # Extract observation time and N
    t_list = list(set(df[key_t]))
    # Time in which scale factor to be introduced
    # When N == 1 (i.e., photometry), no scale factor is introduced.
    t_cor_list = []
    for t in t_list:
        df_t = df[df[key_t] == t]
        N_t = len(df_t)
        #print(f"{key_t}={t} N={N_t}")
        # Do not introduce scale factor when N = 1 (i.e., photometry)
        if N_t > 1:
            t_cor_list.append(t)
        else:
            # Add into list of updated dataframe
            df1_list.append(df_t)

    N_sf = len(t_cor_list)
    print(f"      Number of scale factors = {N_sf}")
     
    # Determine scale factor
    # List of scale factor to be searched
    sf_list = np.arange(sf0, sf1+sfstep, sfstep)
    for idx_t, t_cor in enumerate(t_cor_list):
        df_t = df[df[key_t] == t_cor].copy()

        # Minimize chi2
        for idx_sf, sf in enumerate(sf_list):
            df_t["scalefactor"] = sf
            if idx_sf == 0:
                # Calculate chi-squared
                # These are bench marks
                chi2 = calc_chi2(df_t)
                scalefactor = sf
            else:
                # Calculate chi2
                chi2_sf = calc_chi2(df_t)
                if chi2_sf < chi2:
                    # Update
                    chi2 = chi2_sf
                    scalefactor = sf
                else:
                    pass
            #print(f"chi2, sf = {chi2:.2f}, {scalefactor}")
        #print(f"{idx_t+1}-th scale factor {scalefactor}")
        df_t["scalefactor"] = scalefactor
        # Do not save for the consistency 
        # (N=1 data has no such column)
        df_t = df_t.drop(["diff"], axis=1)

        df1_list.append(df_t)

    # Make a new dataframe with updated scale factors
    df1 = pd.concat(df1_list)
    # Sort by time
    df1 = df1.sort_values(by=[key_t])
    return df1
