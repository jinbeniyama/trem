#!/usr/bin/bash
# Code for TPM test

# Input arguments =============================================================
## shape model
fOBJ=$1
## spin file
fSPIN=$2
## obs file
fOBS=$3
## ephem file
fEPH=$4

# Parameters ==================================================================
## Emissivity
EPS=0.9
## Bond albedo
BondA=0.12
TI=100
## Roughness (crater) parameters
## See Hung+2022, PSJ for the details.
CA=50
CR=0.5
# Parameters ==================================================================

# Run TPM
# wo/-f: Lagerros approximation 
# w/-f : Full heat diffusion within the craters, much slower, but safer, more physical
# Just output the command
#echo ${fOBJ} ${fEPH} ${EPS} ${TI} ${BondA} ${CA} ${CR} | runtpm -f -S ${fSPIN} -o ${fOBS} > ${OUTDIR}/tpmout_ti${TI}_ca${CA}_cr${CR}.dat &
echo ${fOBJ} ${fEPH} ${EPS} ${TI} ${BondA} ${CA} ${CR} | runtpm -S ${fSPIN} -o ${fOBS} > ./tpmout_ti${TI}_ca${CA}_cr${CR}.dat
