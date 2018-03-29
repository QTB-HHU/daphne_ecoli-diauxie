#!/bin/bash 
##################################################################################
# Author: Antonella Succurro
#         a.succurro AT uni-koeln DOT de
###
# Description: reproduce Figure 2 from Succurro et al. 2018
###
# Depends on data in ../../ecoli/enjalbert2015_data/
###
# Usage: --runglucose reproduce conditions 15mM Glc as in Enjalbert et al. 2015 Fig. 2A
#        --runhighacetate  reproduce conditions 45mM Ac
#        --runsingle --runecgl --ratioecgl "1." runs a single E. coli model
#        -M runs MOMA, else pFBA
##################################################################################
OUTDIR=../../outputs/fig2_enjalbert2015_geneexp/
mkdir -p ${OUTDIR}
source ../../pubvenv/bin/activate

maxsteps=10000
ibm=0.0027
rxns="ACKr PPCK PPS FBP ICL MALS PFK PYK PPC ICDHyr"
rxns0="ACKr PPCK PPS FBP ICL MALS PFK PYK"
MOMAPFBADIR=${OUTDIR}/../fig1_enjalbert2015_fig2a/

python ecDiauxie.py -b ${ibm} -n ${maxsteps} -t 1. -x '-11.5' -p ${OUTDIR} --run --runsingle --runecgl --ratioecgl "1." --runglucose -l pfba_notrans_mono
python ecDiauxie.py -b ${ibm} -n ${maxsteps} -t 1. -x '-11.5' -p ${OUTDIR} --run --runsingle --runecgl --ratioecgl "1." --runhighacetate -l pfba_notrans_mono
### Fig. 2a
python analyseFluxes.py -i "${OUTDIR}/endOfSimulation-ecoli_core-batch_high_Ac-pfba_notrans_mono_ECgl.p ${OUTDIR}/endOfSimulation-ecoli_core-batch_low_Glc-pfba_notrans_mono_ECgl.p" -r "${rxns0}" -C -l 'Acetate Glucose' -p ${OUTDIR}
### Fig. 2b
python analyseFluxes.py -i "${MOMAPFBADIR}/endOfSimulation-ecoli_core-batch_low_Glc-moma_notrans_mono_ECgl.p ${MOMAPFBADIR}/endOfSimulation-ecoli_core-batch_low_Glc-pfba_notrans_mono_ECgl.p" -C -l 'MOMA pFBA' -p ${OUTDIR} -r "${rxns}" -G -k 0.3


OUTDIR=../../outputs/supplementary_figures/geneexp/
mkdir -p ${OUTDIR}
### Fig. S5
python analyseFluxes.py -i "${MOMAPFBADIR}/endOfSimulation-ecoli_core-batch_low_Glc-moma_notrans_mono_ECgl.p ${MOMAPFBADIR}/endOfSimulation-ecoli_core-batch_low_Glc-pfba_notrans_mono_ECgl.p" -C -l 'MOMA pFBA' -p ${OUTDIR} -r "${rxns}" -G -k 0.5



deactivate
