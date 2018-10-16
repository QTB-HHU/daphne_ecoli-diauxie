#!/bin/bash
##################################################################################
# Author: Antonella Succurro
#         a.succurro AT uni-koeln DOT de
###
# Description: reproduce Figure 1 from Succurro et al. 2018
###
# Depends on data in ../../ecoli/varma1994_data/
###
# Usage: --runvarmabatch runs conditions as in Varma and Palsson 1994 Fig. 7
#        --runvarmafedbatch runs conditions as in Varma and Palsson 1994 Fig. 10
#        -M runs MOMA, else pFBA
##################################################################################

OUTDIR=../../outputs/supplementary_figures/varma1994/
mkdir -p ${OUTDIR}
source ../../pubvenv/bin/activate

maxsteps=10000

## Fig 7 (batch) 0.01 gDW/L: for vol 30mL  ibm = 0.0003
## Fig 10 (fed) 0.008 gDW/L: for vol 30mL  ibm = 0.00024
#ibm=0.002

### Simulate batch and fed-batch with MOMA and pFBA
ibm=0.0003
python ecDiauxie.py -b ${ibm} -n ${maxsteps} -t 10 -x '-11.5' -p ${OUTDIR} --run --runsingle --runecgl --ratioecgl "1." --runvarmabatch -l pfba_notrans_mono
python ecDiauxie.py -b ${ibm} -n ${maxsteps} -t 10 -x '-11.5' -p ${OUTDIR} --run --runsingle --runecgl --ratioecgl "1." --runvarmabatch -M -l moma_notrans_mono
ibm=0.00024
python ecDiauxie.py -b ${ibm} -n ${maxsteps} -t 10 -x '-11.5' -p ${OUTDIR} --run --runsingle --runecgl --ratioecgl "1." --runvarmafedbatch -l pfba_notrans_mono
python ecDiauxie.py -b ${ibm} -n ${maxsteps} -t 10 -x '-11.5' -p ${OUTDIR} --run --runsingle --runecgl --ratioecgl "1." --runvarmafedbatch -M -l moma_notrans_mono

### Produce plot
python ecDiauxie.py --runsingle --runvarmabatch  -p ${OUTDIR}  -m ${OUTDIR}/endOfSimulation-ecoli_core-varma_batch-pfba_notrans_mono_ECgl.p
python ecDiauxie.py --runsingle --runvarmabatch  -p ${OUTDIR}  -m ${OUTDIR}/endOfSimulation-ecoli_core-varma_batch-moma_notrans_mono_ECgl.p -M
python ecDiauxie.py --runsingle --runvarmafedbatch  -p ${OUTDIR}  -m ${OUTDIR}/endOfSimulation-ecoli_core-varma_fedbatch-pfba_notrans_mono_ECgl.p
python ecDiauxie.py --runsingle --runvarmafedbatch  -p ${OUTDIR}  -m ${OUTDIR}/endOfSimulation-ecoli_core-varma_fedbatch-moma_notrans_mono_ECgl.p -M

python agreement.py -p ${OUTDIR} -m ${OUTDIR}/endOfSimulation-ecoli_core-varma_batch-pfba_notrans_mono_ECgl.p
python agreement.py -p ${OUTDIR} -m ${OUTDIR}/endOfSimulation-ecoli_core-varma_fedbatch-pfba_notrans_mono_ECgl.p
python agreement.py -p ${OUTDIR} -m ${OUTDIR}/endOfSimulation-ecoli_core-varma_batch-moma_notrans_mono_ECgl.p
python agreement.py -p ${OUTDIR} -m ${OUTDIR}/endOfSimulation-ecoli_core-varma_fedbatch-moma_notrans_mono_ECgl.p

# General plots
python analyseFluxes.py -i ${OUTDIR}endOfSimulation-ecoli_core-varma_batch-pfba_notrans_mono_ECgl.p -p ${OUTDIR} -T -D -P
python analyseFluxes.py -i ${OUTDIR}endOfSimulation-ecoli_core-varma_batch-moma_notrans_mono_ECgl.p -p ${OUTDIR} -T -D -P
python analyseFluxes.py -i ${OUTDIR}endOfSimulation-ecoli_core-varma_fedbatch-pfba_notrans_mono_ECgl.p -p ${OUTDIR} -T -D -P
python analyseFluxes.py -i ${OUTDIR}endOfSimulation-ecoli_core-varma_fedbatch-moma_notrans_mono_ECgl.p -p ${OUTDIR} -T -D -P

# Exchange fluxes
python analyseFluxes.py -i ${OUTDIR}endOfSimulation-ecoli_core-varma_batch-pfba_notrans_mono_ECgl.p -p ${OUTDIR} -E
python analyseFluxes.py -i ${OUTDIR}endOfSimulation-ecoli_core-varma_batch-moma_notrans_mono_ECgl.p -p ${OUTDIR} -E
python analyseFluxes.py -i ${OUTDIR}endOfSimulation-ecoli_core-varma_fedbatch-pfba_notrans_mono_ECgl.p -p ${OUTDIR} -E
python analyseFluxes.py -i ${OUTDIR}endOfSimulation-ecoli_core-varma_fedbatch-moma_notrans_mono_ECgl.p -p ${OUTDIR} -E


deactivate
