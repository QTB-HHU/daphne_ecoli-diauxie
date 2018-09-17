#!/bin/bash
##################################################################################
# Author: Antonella Succurro
#         a.succurro AT uni-koeln DOT de
###
# Description: checks for revision
###
# Depends on data in ../../ecoli/varma1994_data/
###
# Usage: --runvarmabatch runs conditions as in Varma and Palsson 1994 Fig. 7
#        --runvarmafedbatch runs conditions as in Varma and Palsson 1994 Fig. 10
#        -M runs MOMA, else pFBA
##################################################################################

OUTDIR=../../outputs/revision_checks_figures/varma1994/
mkdir -p ${OUTDIR}
source ../../pubvenv/bin/activate

maxsteps=10000

## Fig 7 (batch) 0.01 gDW/L: for vol 30mL  ibm = 0.0003
## Fig 10 (fed) 0.008 gDW/L: for vol 30mL  ibm = 0.00024


### PFBA

# ### Simulate batch and fed-batch with pFBA varying EX_O2 uptake and EX_AC secretion
# ### NB!! TO CHECK IN THE CODE THAT -v "${a}" sets upper bound
# for o in 12 11.5 11 10.5 10;
# do
#     sto=${o//\./p}
#     for a in 2 2.5 3 3.5 4 4.5;
#     do
# 	sta=${a//\./p}
# 	ibm=0.0003
# 	python ecDiauxie_revision.py -b ${ibm} -n ${maxsteps} -t 10 -x "-${o}" -v "${a}" -p ${OUTDIR} --run --runsingle --runecgl --ratioecgl "1." --runvarmabatch -l varma_o2_${sto}_ac_${sta}
# 	ibm=0.00024
# 	python ecDiauxie_revision.py -b ${ibm} -n ${maxsteps} -t 10 -x "-${o}" -v "${a}" -p ${OUTDIR} --run --runsingle --runecgl --ratioecgl "1." --runvarmafedbatch -l varma_o2_${sto}_ac_${sta}

# ### Produce plot
# 	python ecDiauxie_revision.py --runsingle --runvarmabatch  -p ${OUTDIR}  -m ${OUTDIR}/endOfSimulation-ecoli_core-varma_batch-varma_o2_${sto}_ac_${sta}_ECgl.p
# 	python ecDiauxie_revision.py --runsingle --runvarmafedbatch  -p ${OUTDIR}  -m ${OUTDIR}/endOfSimulation-ecoli_core-varma_fedbatch-varma_o2_${sto}_ac_${sta}_ECgl.p

#     done
# done


### Simulate batch and fed-batch with pFBA varying EX_ac uptake
### NB!! -q "${a}" sets max lower bound
for a in 3 4 5 6 7 8 9 10;
    do
	sta=${a//\./p}
	ibm=0.0003
	python ecDiauxie_revision.py -b ${ibm} -n ${maxsteps} -t 10 -x "-11.5" -q "${a}" -p ${OUTDIR} --run --runsingle --runecgl --ratioecgl "1." --runvarmabatch -l varma_ac_upt_${sta}
	ibm=0.00024
	python ecDiauxie_revision.py -b ${ibm} -n ${maxsteps} -t 10 -x "-11.5" -q "${a}" -p ${OUTDIR} --run --runsingle --runecgl --ratioecgl "1." --runvarmafedbatch -l varma_ac_upt_${sta}

### Produce plot
	python ecDiauxie_revision.py --runsingle --runvarmabatch  -p ${OUTDIR}  -m ${OUTDIR}/endOfSimulation-ecoli_core-varma_batch-varma_ac_upt_${sta}_ECgl.p
	python ecDiauxie_revision.py --runsingle --runvarmafedbatch  -p ${OUTDIR}  -m ${OUTDIR}/endOfSimulation-ecoli_core-varma_fedbatch-varma_ac_upt_${sta}_ECgl.p
done



##### MOMA

#BATCH
#ibm=0.0003
#python ecDiauxie_revision.py -b ${ibm} -n ${maxsteps} -t 10 -x '-11.5' -p ${OUTDIR} --run --runsingle --runecgl --ratioecgl "1." --runvarmabatch -M -l moma_notrans_mono

#FED-BATCH
#ibm=0.00024
#python ecDiauxie_revision.py -b ${ibm} -n ${maxsteps} -t 10 -x '-11.5' -p ${OUTDIR} --run --runsingle --runecgl --ratioecgl "1." --runvarmafedbatch -M -l moma_notrans_mono

#python ecDiauxie_revision.py --runsingle --runvarmabatch  -p ${OUTDIR}  -m ${OUTDIR}/endOfSimulation-ecoli_core-varma_batch-moma_notrans_mono_ECgl.p -M

#python ecDiauxie_revision.py --runsingle --runvarmafedbatch  -p ${OUTDIR}  -m ${OUTDIR}/endOfSimulation-ecoli_core-varma_fedbatch-moma_notrans_mono_ECgl.p -M

deactivate
