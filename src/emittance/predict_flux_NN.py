#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Predict fluxes using NN model.
"""
import os 
from argparse import ArgumentParser as ap
import numpy as np
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
        "--TI0", type=float, default=10,
        help="Minimum TI")
    parser.add_argument(
        "--TI1", type=float, default=2500,
        help="Maximum TI")
    parser.add_argument(
        "--TIstep", type=float, default=5,
        help="Step pf TI")
    parser.add_argument(
        "--theta0", type=float, default=0,
        help="Minimum theta bar")
    parser.add_argument(
        "--theta1", type=float, default=60,
        help="Maximum theta bar")
    parser.add_argument(
        "--thetastep", type=float, default=1,
        help="Step pf theta bar")
    parser.add_argument(
        "--notuse", type=float, nargs="*", default=None,
        help="Epoch not used")
    parser.add_argument(
        "--outdir", type=str, default="NNprediction",
        help="output file name")
    args = parser.parse_args()

    outdir = args.outdir
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

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

    TI_list = np.arange(args.TI0, args.TI1+args.TIstep, args.TIstep)
    theta_list = np.arange(args.theta0, args.theta1+args.thetastep, args.thetastep)
    
    for i in range(len(epoch_unique_array)):
    
        # Get list of wavelengths (wave_list)
        dist = np.abs(TPM_sims[:,5]-epoch_unique_array[i])
        idx1 = np.where(dist<1e-8)
        TPM_sims_epoch = TPM_sims[idx1[0],:]
        idx2 = np.where(np.logical_and(TPM_sims_epoch[:,0]==TI0,TPM_sims_epoch[:,1]==theta0))[0]
        wave_list = TPM_sims_epoch[idx2,6]
    
        # Get NN prediction for i-th epoch
        # Length is N (=AxBxC), where A is len(wave_list), B is len(TI_list), and C is len(theta_list)
        A = len(wave_list)
        B = len(TI_list)
        C = len(theta_list)
        N = A * B * C
        print(f"For epoch {i+1}, N = A x B x C = {A} x {B} x {C} = {N}")
    
        W, TI, THETA = np.meshgrid(wave_list, TI_list, theta_list, indexing='ij')
    
        wave_list4query = W.ravel()
        TI_list4query = TI.ravel()
        theta_list4query = THETA.ravel()
    
        f_model = predict_flux_NN(
            args.modeldir, epoch_unique_array[i], wave_list4query, 
            TI_list4query, theta_list4query)


        shape = f_model.shape
        epoch_arr = np.full(shape, epoch_unique_array[i], dtype=float)
        arr = np.stack([
            epoch_arr,
            TI_list4query,
            theta_list4query,
            wave_list4query,
            f_model,
        ], axis=0)

        outdir = args.outdir
        # Save as "LUT_2450991.767627034.npy"
        f_out = f"LUT_{epoch_unique_array[i]}.npy"
        out = os.path.join(outdir, f_out)
        np.save(out, arr.T)
