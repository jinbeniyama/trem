#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Make look-up-table to make NN model.
"""
import os 
from argparse import ArgumentParser as ap
import pandas as pd
import numpy as np
import glob
import re


def parseFile(fname, outfile):
    with open(fname, 'r') as file:
        data = file.read()

    ### extract TPM simualtion parameters
    r0=data.find('r>')+2
    r1=data.find('\n',r0)
    paramList=re.findall(r"[-+]?(?:\d*\.*\d+)", data[r0:r1])
    #print(paramList)
    gamma = paramList[0]
    tBar = paramList[1]
    gammaC = paramList[2]
    gammaRho = paramList[3]
    bondA = paramList[4]

    ### extract TPM simualtion fluxes
    r0=0
    rMax = len(data)

    while (r0<rMax):
        r0=data.find('f>',r0)
        r1=data.find('\n',r0)
        #print(r0, r1, rMax)
        if ((r1<0) | (r0<0)):
            break;
        if ((r0>0) & (r1>0)):
            paramList=re.findall(r"[-+]?(?:\d*\.*\d+)", data[r0:r1])
        #print(paramList)
            jd =  paramList[1]
            wave = paramList[3]
            tpmFlux = paramList[6]
            print(gamma, tBar, gammaC, gammaRho, bondA, jd, wave, tpmFlux, sep=",", file=outfile)
            r0=r1+1


if __name__ == "__main__":
    parser = ap(
        description="Make look-up-table.")
    parser.add_argument(
        "res", type=str,
        help="TPMresult")
    parser.add_argument(
        "--out", type=str, default="lut.txt",
        help="output file name")
    args = parser.parse_args()

    
    resall = glob.glob(f'{args.res}/tpmout*')
    N_res = len(resall)
    print(f"N_res = {N_res}")
    assert N_res == 440, "Check TPM results."
    
    outf=open(args.out, "w")
    for fn in resall:
        try:
            parseFile(fn,outf)
        except:
            print("error on ", fn)
    outf.close()
