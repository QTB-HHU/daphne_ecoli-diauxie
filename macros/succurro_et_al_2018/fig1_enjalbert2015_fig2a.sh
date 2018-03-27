#!/bin/bash 
##################################################################################
# Author: Antonella Succurro
#         a.succurro AT uni-koeln DOT de
###
# Description: reproduce Figure 1 from Succurro et al. 2018
###
# Depends on data in ../../ecoli/enjalbert2015_data/
###
# Usage: --runglucose reproduce conditions 15mM Glc as in Enjalbert et al. 2015 Fig. 2A
#        --runfedlowacetate  reproduce conditions 15mM Glc & constant 4mM Ac  as in Enjalbert et al. 2015 Fig. 6A
#        --runsingle --runecgl --ratioecgl "1." runs a single E. coli model
#        -M runs MOMA, else pFBA
##################################################################################
OUTDIR=../../outputs/fig1_enjalbert2015_fig2a/
mkdir -p ${OUTDIR}
source ../../pubvenv/bin/activate

maxsteps=10000
ibm=0.0026

### Simulate 15mM Glc condition with MOMA and pFBA
python ecDiauxie.py -b ${ibm} -n ${maxsteps} -t 10. -x '-11.5' -p ${OUTDIR} --run --runsingle --runecgl --ratioecgl "1." --runglucose -l pfba_notrans_mono
python ecDiauxie.py -b ${ibm} -n ${maxsteps} -t 10 -x '-11.5' -p ${OUTDIR} --run --runsingle --runecgl --ratioecgl "1." --runglucose -M -l moma_notrans_mono

### Produce plots
python ecDiauxie.py --runsingle --runecgl  -p ${OUTDIR}  -m ${OUTDIR}/endOfSimulation-ecoli_core-batch_low_Glc-pfba_notrans_mono_ECgl.p
python ecDiauxie.py --runsingle --runecgl  -p ${OUTDIR}  -m ${OUTDIR}/endOfSimulation-ecoli_core-batch_low_Glc-moma_notrans_mono_ECgl.p -M

python agreement.py -p ${OUTDIR} -m ${OUTDIR}/endOfSimulation-ecoli_core-batch_low_Glc-pfba_notrans_mono_ECgl.p
python agreement.py -p ${OUTDIR} -m ${OUTDIR}/endOfSimulation-ecoli_core-batch_low_Glc-moma_notrans_mono_ECgl.p


deactivate
