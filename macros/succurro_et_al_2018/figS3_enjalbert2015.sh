#!/bin/bash
OUTDIR=../../outputs/supplementary_figures/enjMOMApFBA/
INPDIR=../../outputs/fig1_enjalbert2015_fig2a/
mkdir -p ${OUTDIR}
conda activate daphnePy2

python analyseFluxes.py -i ${INPDIR}endOfSimulation-ecoli_core-batch_low_Glc-pfba_notrans_mono_ECgl.p -p ${OUTDIR} -T -D -P
python analyseFluxes.py -i ${INPDIR}endOfSimulation-ecoli_core-batch_low_Glc-moma_notrans_mono_ECgl.p -p ${OUTDIR} -T -D -P


