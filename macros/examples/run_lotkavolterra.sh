#!/bin/bash
OUTDIR=../../outputs/lotkavolterra/
mkdir -p ${OUTDIR}
source ../../pubvenv/bin/activate
python lotkavolterra_pureODE.py -p ${OUTDIR}
deactivate
