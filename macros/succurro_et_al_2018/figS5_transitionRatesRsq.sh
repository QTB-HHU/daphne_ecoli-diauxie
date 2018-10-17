#!/bin/bash 
##################################################################################
# Author: Antonella Succurro
#         a.succurro AT uni-koeln DOT de
###
# Description: reproduce Figure 5b from Succurro et al. 2018
###
# Depends on data in ../../ecoli/enjalbert2015_data/
###
# Usage: --runfedhighacetate  reproduce conditions 15mM Glc & constant ~32mM Ac  as in Enjalbert et al. 2015 Fig. 6C
#        --runfedlowacetate  reproduce conditions 15mM Glc & constant 4mM Ac  as in Enjalbert et al. 2015 Fig. 6A
#        --runconsortium --ratioecgl "0.75" runs two E. coli models with starting biomass ratio for ECgl =  0.75 TOT
#        -e '0.9' --phitransition --psitransition activates ECgl <-> ECac transitions, with efficiency 0.9
##################################################################################
OUTDIR=../../outputs/supplementary_figures/
mkdir -p ${OUTDIR}
source ../../pubvenv/bin/activate

maxsteps=10000
minstep=0.0

##################### Constant rates
psio=0.04
phio=0.04
vpsi=0.0
vphi=0.0
kpsi=0.0
kphi=0.0

OUTDIR=../../outputs/supplementary_figures/constant_rates_figs4-5a-5b/
mkdir -p ${OUTDIR}

ibm=0.0027
python ecDiauxie.py -b ${ibm} -n ${maxsteps} -t 10 -x '-11.5' -s ${minstep} -p ${OUTDIR}  --run --runconsortium --ratioecgl "0.95" --runglucose -l pfba_wphipsitrans_two -e '0.9' --phitransition --psitransition --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}"
python ecDiauxie.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-batch_low_Glc-pfba_wphipsitrans_two_ECgl_ECac_0p95.p --psitransition --phitransition --ratioecgl "0.95"

ibm=0.006
python ecDiauxie.py -b ${ibm} -n ${maxsteps} -t 11 -x '-11.5' -s ${minstep} -p ${OUTDIR} --run --runconsortium --ratioecgl "0.75" --runfedhighacetate -l pfba_wphipsitrans_two -e '0.9' --phitransition --psitransition --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}"
python ecDiauxie.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-fedbatch_high_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p75.p --phitransition --psitransition --ratioecgl "0.75"

ibm=0.0038
python ecDiauxie.py -b ${ibm} -n ${maxsteps} -t 13 -x '-11.5' -s ${minstep} -p ${OUTDIR} --run --runconsortium --ratioecgl "0.95" --runfedlowacetate -l pfba_wphipsitrans_two -e '0.9' --phitransition --psitransition  --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}"
python ecDiauxie.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-fedbatch_low_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p95.p --ratioecgl "0.95" --phitransition --psitransition 
#python ecDiauxie.py -b ${ibm} -n ${maxsteps} -t 13 -x '-11.5' -s ${minstep} -p ${OUTDIR} --run --runconsortium --ratioecgl "0.75" --runfedlowacetate -l pfba_wphipsitrans_two -e '0.9' --phitransition --psitransition  --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}"
#python ecDiauxie.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-fedbatch_low_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p75.p --ratioecgl "0.75" --phitransition --psitransition 

######################

psio=0.04
phio=0.04
kpsi=30.0
kphi=5.0
vpsi=0.2
vphi=0.2

OUTDIR=../../outputs/supplementary_figures/mmrates_figs4-5a-5b/
mkdir -p ${OUTDIR}

ibm=0.0027
python ecDiauxie.py -b ${ibm} -n ${maxsteps} -t 10 -x '-11.5' -s ${minstep} -p ${OUTDIR}  --run --runconsortium --ratioecgl "0.95" --runglucose -l pfba_wphipsitrans_two -e '0.9' --phitransition --psitransition --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}"
python ecDiauxie.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-batch_low_Glc-pfba_wphipsitrans_two_ECgl_ECac_0p95.p --psitransition --phitransition --ratioecgl "0.95"

ibm=0.006
python ecDiauxie.py -b ${ibm} -n ${maxsteps} -t 11 -x '-11.5' -s ${minstep} -p ${OUTDIR} --run --runconsortium --ratioecgl "0.75" --runfedhighacetate -l pfba_wphipsitrans_two -e '0.9' --phitransition --psitransition --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}"
python ecDiauxie.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-fedbatch_high_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p75.p --phitransition --psitransition --ratioecgl "0.75"

ibm=0.0038
python ecDiauxie.py -b ${ibm} -n ${maxsteps} -t 13 -x '-11.5' -s ${minstep} -p ${OUTDIR} --run --runconsortium --ratioecgl "0.95" --runfedlowacetate -l pfba_wphipsitrans_two -e '0.9' --phitransition --psitransition  --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}"
python ecDiauxie.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-fedbatch_low_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p95.p --ratioecgl "0.95" --phitransition --psitransition 
# python ecDiauxie.py -b ${ibm} -n ${maxsteps} -t 13 -x '-11.5' -s ${minstep} -p ${OUTDIR} --run --runconsortium --ratioecgl "0.75" --runfedlowacetate -l pfba_wphipsitrans_two -e '0.9' --phitransition --psitransition  --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}"
# python ecDiauxie.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-fedbatch_low_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p75.p --ratioecgl "0.75" --phitransition --psitransition 

######################
OUTDIR=../../outputs/supplementary_figures/mmrates_figs4-5a-5b/
echo "with MM rates"
python agreement.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-batch_low_Glc-pfba_wphipsitrans_two_ECgl_ECac_0p95.p
python agreement.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-fedbatch_high_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p75.p
python agreement.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-fedbatch_low_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p95.p
#python agreement.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-fedbatch_low_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p75.p

######################
OUTDIR=../../outputs/supplementary_figures/constant_rates_figs4-5a-5b/
echo "with constant rates"
python agreement.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-batch_low_Glc-pfba_wphipsitrans_two_ECgl_ECac_0p95.p
python agreement.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-fedbatch_high_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p75.p
python agreement.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-fedbatch_low_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p95.p
#python agreement.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-fedbatch_low_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p75.p

deactivate
