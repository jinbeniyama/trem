#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Predict fluxes using NN model.
"""
import os 
from argparse import ArgumentParser as ap
import numpy as np
import pandas as pd
import keras
import pickle


def predict_flux_NN(modeldir, Epoch, wavelength, Gamma, theta):
    # Function to evaluate neural networks saved in "saved_model" folder
    #
    # Inputs are:
    # - Epoch: one of the unique 22 epochs in LUT20250508.txt, float
    # - Wavelength: Nx1 array of wavelengths in microns
    # - Gamma: Nx1 array of thermal inertia in J m-2 s-1/2 k-1
    # - Theta: Nx1 array of hapke angle in degrees
    #
    # Outputs are:
    #
    # - Flux: Nx1 array of flux predictions in Jy
    #
    # Note 1: Wavelength, Gamma and Theta must be between their min and max values in LUT20250508
    # Note 2: Wavelength, Gamma and Theta must be single-value floats or array of the same length

    # Load the model 
    # before 2025-09-13
    #model = keras.models.load_model(modeldir+'/NN_'+str(Epoch)+'.keras')
    #with open(modeldir+r"/scaler_in_"+str(Epoch)+".pkl", "rb") as input_file:
    #    scaler = pickle.load(input_file)

    # after 2025-09-13 (just the filenames are updated)
    model = keras.models.load_model(modeldir+'/NN_'+str(Epoch))
    with open(modeldir+r"/scaler_"+str(Epoch)+".pkl", "rb") as input_file:
        scaler = pickle.load(input_file)

    # Prepare the input array

    if isinstance(wavelength, float) or isinstance(wavelength, int): # inputs are not arrays

        X_input = np.vstack([[Gamma,theta,wavelength],[Gamma,theta,wavelength]])

        # Predict the flux
        Flux = model.predict(scaler.transform(X_input),verbose = 0)[0][0].flatten()

    else: # inputs are all arrays

        X_input = np.vstack([[Gamma,theta,wavelength]]).T
        # Predict the flux
        Flux = model.predict(scaler.transform(X_input),verbose = 0).flatten()

    return Flux


def add_noise_by_epoch(df_model, df_obs):
    """
    jd ごとにモデルと観測を対応させてノイズを加える
    df_model: columns ['jd', 'w', 'f_obs']
    df_obs: columns ['jd', 'wavelength', 'fluxerr']
    """
    result_list = []

    # ユニークな jd をループ
    for jd_epoch in df_model['jd'].unique():
        df_mod_jd = df_model[df_model['jd'] == jd_epoch].copy()
        df_obs_jd = df_obs[df_obs['jd'] == jd_epoch].copy()

        if df_obs_jd.empty:
            # 観測がない場合は NaN
            df_mod_jd['ferr_obs'] = np.nan
            df_mod_jd['f_obs'] = np.nan
            result_list.append(df_mod_jd[['jd','w','f_obs','ferr_obs']])
            continue

        # 波長で sort
        df_mod_jd = df_mod_jd.sort_values('w').reset_index(drop=True)
        df_obs_jd = df_obs_jd.sort_values('wavelength').reset_index(drop=True)

        # 観測数とモデル数が一致しない場合は min(n) で対応
        n_match = min(len(df_mod_jd), len(df_obs_jd))
        df_mod_jd = df_mod_jd.iloc[:n_match]
        df_obs_jd = df_obs_jd.iloc[:n_match]

        # ノイズ追加
        noise = np.random.normal(loc=0, scale=df_obs_jd['fluxerr'].values)
        df_mod_jd['f_obs'] = df_mod_jd['f_obs'].values + noise
        df_mod_jd['ferr_obs'] = df_obs_jd['fluxerr'].values

        result_list.append(df_mod_jd[['jd','w','f_obs','ferr_obs']])

    # 全て結合
    df_tmp = pd.concat(result_list, ignore_index=True)
    return df_tmp

if __name__ == "__main__":
    parser = ap(
        description="Predict fluxes using NN model.")
    parser.add_argument(
        "lut", type=str,
        help="Look-up-table make with 'make_lut.py'")
    parser.add_argument(
        "modeldir", type=str,
        help="Directory with NN model")
    parser.add_argument(
        "obsflux", type=str,
        help="Observations (to extract error)")
    parser.add_argument(
        "--TIrego", type=float, default=50,
        help="Minimum TI")
    parser.add_argument(
        "--TIrock", type=float, default=1000,
        help="Maximum TI")
    parser.add_argument(
        "--Htheta", type=float, default=30,
        help="Theta bar")
    parser.add_argument(
        "--alpha", type=float, default=0.5,
        help="Regolith abundance")
    parser.add_argument(
        "--notuse", type=float, nargs="*", default=None,
        help="Epoch not used")
    parser.add_argument(
        "--wonoise", action="store_true", default=False,
        help="Do not add noise to simulation data")
    parser.add_argument(
        "--out", type=str, default="simu_flux.txt",
        help="output file name")
    args = parser.parse_args()
    
    # Maybe useless. Either lut or obsflux is fine.
    TPM_sims = np.loadtxt(args.lut, delimiter = ',')
    epoch_unique_array = np.unique(TPM_sims[:,5])

    print(f"Unique epochs: N={len(epoch_unique_array)}")

    # Remove useless epochs here
    if args.notuse:
        for epoch_notuse in args.notuse:
            N0 = len(epoch_unique_array) 
            epoch_unique_array = [x for x in epoch_unique_array if x != epoch_notuse]
            N1 = len(epoch_unique_array) 
            if N1-N0 != 0:
                print(f"Epoch {epoch_notuse} is removed.")

        print(f"  Updated unique epochs: N={len(epoch_unique_array)}")
    
    # A thermal inertia and a Hapke thetabar.
    # These are just to extract wavelength
    line0  = TPM_sims[0]
    TI0    = line0[0]
    theta0 = line0[1]

    df_obs = pd.read_csv(args.obsflux, sep=" ")
    TIrego = args.TIrego
    TIrock = args.TIrock
    Htheta = args.Htheta
    alpha = args.alpha

    df_list = []
    for i in range(len(epoch_unique_array)):
    
        # Get list of wavelengths (wave_list)
        dist = np.abs(TPM_sims[:,5]-epoch_unique_array[i])
        idx1 = np.where(dist<1e-8)
        TPM_sims_epoch = TPM_sims[idx1[0],:]
        idx2 = np.where(np.logical_and(TPM_sims_epoch[:,0]==TI0,TPM_sims_epoch[:,1]==theta0))[0]
        wave_list = TPM_sims_epoch[idx2,6]
    
        # Get NN prediction for i-th epoch
        A = len(wave_list)
        print(f"For epoch {i+1}, N_wave = {A}")
        
        # Regolith component
        W, TI, THETA = np.meshgrid(wave_list, [TIrego], [Htheta], indexing='ij')
        wave_list4query = W.ravel()
        TI_list4query = TI.ravel()
        theta_list4query = THETA.ravel()
        f_model_rego = predict_flux_NN(
            args.modeldir, epoch_unique_array[i], wave_list4query, 
            TI_list4query, theta_list4query)

        # Rock component
        W, TI, THETA = np.meshgrid(wave_list, [TIrock], [Htheta], indexing='ij')
        wave_list4query = W.ravel()
        TI_list4query = TI.ravel()
        theta_list4query = THETA.ravel()
        f_model_rock = predict_flux_NN(
            args.modeldir, epoch_unique_array[i], wave_list4query, 
            TI_list4query, theta_list4query)

        # Sum of the model
        f_model_sum = alpha*f_model_rego + (1-alpha)*f_model_rock

        # Make a dataframe
        # Save f_model as f_obs!
        df_e = pd.DataFrame({
            "w": wave_list,
            "f_obs": f_model_sum
            })
        df_e["jd"] = epoch_unique_array[i]


        df_list.append(df_e)

    df = pd.concat(df_list)

    assert not df.isna().any().any(), "NaN values detected!"

    # Add noise
    print()
    if args.wonoise:
        print("Do not add noise!")
        # Dummy
        df["ferr_obs"] = 0.1
    else:
        print("Add noise using fobs_err.")
        df = add_noise_by_epoch(df, df_obs)

        df_nan = df[df.isna().any(axis=1)]
        assert not df.isna().any().any(), "NaN values detecded! (after adding noise)"

    out = args.out
    df.to_csv(out, sep=" ", index=False)
