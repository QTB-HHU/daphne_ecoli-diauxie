#!/bin/bash 
##################################################################################
# Author: Antonella Succurro
#         a.succurro AT uni-koeln DOT de
###
# Description: reproduce Figure 4 from Succurro et al. 2018
###
# Depends on data in ../../ecoli/enjalbert2015_data/
###
# Usage: --runfedhighacetate  reproduce conditions 15mM Glc & constant ~32mM Ac  as in Enjalbert et al. 2015 Fig. 6C
#        --runfedlowacetate  reproduce conditions 15mM Glc & constant 4mM Ac  as in Enjalbert et al. 2015 Fig. 6A
#        --runconsortium --ratioecgl "0.75" runs two E. coli models with starting biomass ratio for ECgl =  0.75 TOT
#        -e '0.9' --phitransition --psitransition activates ECgl <-> ECac transitions, with efficiency 0.9
##################################################################################
OUTDIR=../../outputs/fig4_enjalbert2015_fig6a_fig6c/
mkdir -p ${OUTDIR}
source ../../pubvenv/bin/activate

maxsteps=10000
ibm=0.0036
minstep=0.0

psio=0.04
phio=0.04
kpsi=30.0
kphi=5.0
vpsi=0.2
vphi=0.2

python ecDiauxie.py -b ${ibm} -n ${maxsteps} -t 13 -x '-11.5' -s ${minstep} -p ${OUTDIR} --run --runconsortium --ratioecgl "1.0" --runfedlowacetate -l pfba_wphipsitrans_two -e '0.9' --phitransition --psitransition --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}"
python ecDiauxie.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-fedbatch_low_Ac-pfba_wphipsitrans_two_ECgl_ECac_1p0.p --ratioecgl "1.0" --phitransition --psitransition 

ibm=0.0058

python ecDiauxie.py -b ${ibm} -n ${maxsteps} -t 11 -x '-11.5' -s ${minstep} -p ${OUTDIR} --run --runconsortium --ratioecgl "0.75" --runfedhighacetate -l pfba_wphipsitrans_two -e '0.9' --phitransition --psitransition --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}"
python ecDiauxie.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-fedbatch_high_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p75.p --phitransition --psitransition --ratioecgl "0.75"

deactivate
