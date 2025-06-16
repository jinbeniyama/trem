#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Merge observations and NN model fluxes.
"""
import os 
from argparse import ArgumentParser as ap
import pandas as pd
import numpy as np


def read_obs(f):
    """Open TPM result and extract observations.

    Parameter
    ---------
    f : str
        output of runtpm

    Return
    ------
    df : pandas.DataFrame
        dataframe with jd, w, f_obs, ferr_obs
    """
    jd_list, w_list = [], []
    f_obs_list, ferr_obs_list = [], []
    with open (f, "r") as f:
        f = f.readlines()
        for l in f:
            if l[0:2] == "f>":
                # Extract info
                l = l.split()
                epoch, jd, n       = l[1], l[2], l[3]
                w, f_obs, ferr_obs = l[4], l[5], l[6]
                jd_list.append(float(jd))
                w_list.append(float(w))
                f_obs_list.append(float(f_obs))
                ferr_obs_list.append(float(ferr_obs))

    df = pd.DataFrame({
        "jd": jd_list,
        "w": w_list,
        "f_obs": f_obs_list,
        "ferr_obs": ferr_obs_list,
        })
    return df


def prepro_NN(f_list):
    """ Read 22 files (JD, TI, Htheta, Wavelength, Predicted flux) for each epoch.
    Make a single DataFrame with JD, TI, Htheta, Wavelength, f_model

    Parameter
    ---------
    f_list : array-like
        list of input files

    Return
    ------
    df : pandas.DataFrame
        merged dataframe
    """
    df_list = []
    for f in f_list:

        # JD, TI, Htheta, Wavelength, Predicted flux
        arr = np.load(f)
        df = pd.DataFrame(arr, columns=["jd", "TI", "Htheta", "w", "f_model"])
        #print(f"N={len(df)}")

        # Round flux
        df["f_model"] = np.round(df["f_model"], 5)
        # To have a consistency
        #df["w"] = np.round(df["w"], 2)

        df_list.append(df)
    df = pd.concat(df_list)
    return df


def check_epoch_wavelength(df, key_epoch="jd", key_w="w"):
    """Check consistency.

    Parameters
    ----------
    df : pandas.DataFrame
        input dataframe
    key_epoch : str
        keyrword for epoch
    key_w : str
        keyword for wavelength

    Returns
    --------
    epoch_unique : array-like
        Unique epoch
    w_unique : array-like
        Unique wavelength
    """
    epoch_unique = list(set(df[key_epoch]))
    epoch_unique = sorted(epoch_unique)

    w_unique = list(set(df[key_w]))
    w_unique = sorted(w_unique)

    return epoch_unique, w_unique


def concat_obs_NN(df_obs, df_NN):
    """Concat NN and obs and Make a dataframe with jd, w, f_obs, ferr_obs, f_model, scalefactor, TI, Htheta, A

    Parameters
    ----------
    df_obs : pandas.DataFrame
        dataframe with observations
    df_NN : pandas.DataFrame
        dataframe with NN predicted fluxes

    Return
    ------
    df_result : pandas.DataFrame
        resultant dataframe
    """
    # use useful columns
    df_obs_reduced = df_obs[["jd", "w", "f_obs", "ferr_obs"]]
    # concat with keys of "jd" and "w"
    df_result = df_NN.merge(df_obs_reduced, on=["jd", "w"], how="left")
    return df_result


if __name__ == "__main__":
    parser = ap(
        description="Merge observations and NN model fluxes.")
    parser.add_argument(
        "NNdir", type=str,
        help="Directory with NN predictions (.npy etc.)")
    parser.add_argument(
        "f_obs", type=str,
        help="Example of TPM result (to extract observations)")
    parser.add_argument(
        "--out", type=str, default="NN_TPMres.txt.txt",
        help="output file name")
    args = parser.parse_args()

    
    # Read NN predictions
    ## directory where the files are located
    ## Check this directory carefully!
    f_list = os.listdir(args.NNdir)
    Nepoch = len(f_list)
    print(f"  Number of epochs: {Nepoch}")
    ## Add directory
    f_list = [f"{args.NNdir}/{x}" for x in f_list]
    
    ## Concat files
    df_NN = prepro_NN(f_list)
    
    ## Check unique epochs and wavelengths
    epoch_NN_unique, w_NN_unique = check_epoch_wavelength(df_NN)
    print(f"  N={len(epoch_NN_unique)} {epoch_NN_unique}")
    print(f"  N={len(w_NN_unique)} {w_NN_unique}")
    
    # Read observations
    ## read a result of TPM
    df_obs = read_obs(args.f_obs)
    
    ## Check epoch and wavelength
    epoch_obs_unique, w_obs_unique = check_epoch_wavelength(df_obs)
    print(f"N={len(epoch_obs_unique)} {epoch_obs_unique}")
    print(f"N={len(w_obs_unique)} {w_obs_unique}")
    
    # Check if the columns are identical
    # This should be empty
    diff1 = [x for x in w_obs_unique if not x in w_NN_unique]
    assert diff1 == [], "Check the code and inputs."
    # This should be empty
    diff2 = [x for x in w_NN_unique if not x in w_obs_unique]
    assert diff2 == [], "Check the code and inputs."
    
    # Concatenate and save
    df_concat = concat_obs_NN(df_obs, df_NN)
    # This dataframe is equivallent to the return of extract_flux in trem/common.py
    df_concat.to_csv(args.out, sep=" ", index=False)
