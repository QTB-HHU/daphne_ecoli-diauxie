#!/bin/bash 
##################################################################################
# Author: Antonella Succurro
#         a.succurro AT uni-koeln DOT de
###
# Description: reproduce Figure 3 from Succurro et al. 2018
###
# Depends on data in ../../ecoli/enjalbert2015_data/
###
# Usage: --runglucose reproduce conditions 15mM Glc as in Enjalbert et al. 2015 Fig. 2A
#        --runfedlowacetate  reproduce conditions 15mM Glc & constant 4mM Ac  as in Enjalbert et al. 2015 Fig. 6A
#        --runconsortium --ratioecgl "1." runs two E. coli models with starting biomass ratio for ECgl =  0.999 TOT
#        -e '0.9' --phitransition --psitransition activates ECgl <-> ECac transitions, with efficiency 0.9
##################################################################################
OUTDIR=../../outputs/presentations_figures/
mkdir -p ${OUTDIR}
source ../../pubvenv/bin/activate

steps=10000
minstep=0.0


## Simulate uniform

### 15mM Glc condition
ibm=0.0027
python ecDiauxie_isme.py -b ${ibm} -n ${steps} -t 10 -x '-11.5' -s ${minstep} -p ${OUTDIR} --run --runsingle --runecgl --ratioecgl "0.95" --runglucose -l pfba_notrans_mono
python ecDiauxie_isme.py --runsingle --runecgl  -p ${OUTDIR}  -m ${OUTDIR}/endOfSimulation-ecoli_core-batch_low_Glc-pfba_notrans_mono_ECgl.p

### Low Ace
ibm=0.0038
python ecDiauxie_isme.py -b ${ibm} -n ${steps} -t 13 -x '-11.5' -s ${minstep} -p ${OUTDIR} --run --runsingle --runecgl --ratioecgl "0.95" --runfedlowacetate -l pfba_notrans_mono
python ecDiauxie_isme.py --runsingle --runecgl -p ${OUTDIR} -m ${OUTDIR}/endOfSimulation-ecoli_core-fedbatch_low_Ac-pfba_notrans_mono_ECgl.p


### High Ace
ibm=0.006
python ecDiauxie_isme.py -b ${ibm} -n ${steps} -t 11 -x '-11.5' -s ${minstep} -p ${OUTDIR} --run --runsingle --runecgl --ratioecgl "0.75" --runfedhighacetate -l pfba_notrans_mono
python ecDiauxie_isme.py --runsingle --runecgl -p ${OUTDIR} -m ${OUTDIR}/endOfSimulation-ecoli_core-fedbatch_high_Ac-pfba_notrans_mono_ECgl.p


## Simulate stoch

psio=0.04
phio=0.04
vpsi=0.0
kpsi=0.0
vphi=0.0
kphi=0.0


### 15mM Glc condition
ibm=0.0027
python ecDiauxie_isme.py -b ${ibm} -n ${steps} -t 10 -x '-11.5' -s ${minstep} -p ${OUTDIR} --run --runconsortium --ratioecgl "0.95" --runglucose -l pfba_wstochtrans_two -e '0.9' --phitransition --psitransition --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}"
python ecDiauxie_isme.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-batch_low_Glc-pfba_wstochtrans_two_ECgl_ECac_0p95.p --psitransition --phitransition --ratioecgl "0.95"

### Low Ace
ibm=0.0038
python ecDiauxie_isme.py -b ${ibm} -n ${steps} -t 13 -x '-11.5' -s ${minstep} -p ${OUTDIR} --run --runconsortium --ratioecgl "0.95" --runfedlowacetate -l pfba_wstochtrans_two -e '0.9' --phitransition --psitransition --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}"
python ecDiauxie_isme.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-fedbatch_low_Ac-pfba_wstochtrans_two_ECgl_ECac_0p95.p --ratioecgl "0.95" --phitransition --psitransition 


### High Ace
ibm=0.006
python ecDiauxie_isme.py -b ${ibm} -n ${steps} -t 11 -x '-11.5' -s ${minstep} -p ${OUTDIR} --run --runconsortium --ratioecgl "0.75" --runfedhighacetate -l pfba_wstochtrans_two -e '0.9' --phitransition --psitransition --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}"
python ecDiauxie_isme.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-fedbatch_high_Ac-pfba_wstochtrans_two_ECgl_ECac_0p75.p --phitransition --psitransition --ratioecgl "0.75"




## Simulate resp trans

psio=0.04
phio=0.04
kpsi=30.0
kphi=5.0
vpsi=0.2
vphi=0.2

### 15mM Glc condition
ibm=0.0027
python ecDiauxie_isme.py -b ${ibm} -n ${steps} -t 10 -x '-11.5' -s ${minstep} -p ${OUTDIR} --run --runconsortium --ratioecgl "0.95" --runglucose -l pfba_wphipsitrans_two -e '0.9' --phitransition --psitransition --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}"
python ecDiauxie_isme.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-batch_low_Glc-pfba_wphipsitrans_two_ECgl_ECac_0p95.p --psitransition --phitransition --ratioecgl "0.95"

### Low Ace
ibm=0.0038
python ecDiauxie_isme.py -b ${ibm} -n ${steps} -t 13 -x '-11.5' -s ${minstep} -p ${OUTDIR} --run --runconsortium --ratioecgl "0.95" --runfedlowacetate -l pfba_wphipsitrans_two -e '0.9' --phitransition --psitransition --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}"
python ecDiauxie_isme.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-fedbatch_low_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p95.p --ratioecgl "0.95" --phitransition --psitransition 


### High Ace
ibm=0.006
python ecDiauxie_isme.py -b ${ibm} -n ${steps} -t 11 -x '-11.5' -s ${minstep} -p ${OUTDIR} --run --runconsortium --ratioecgl "0.75" --runfedhighacetate -l pfba_wphipsitrans_two -e '0.9' --phitransition --psitransition --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}"
python ecDiauxie_isme.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-fedbatch_high_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p75.p --phitransition --psitransition --ratioecgl "0.75"

deactivate
