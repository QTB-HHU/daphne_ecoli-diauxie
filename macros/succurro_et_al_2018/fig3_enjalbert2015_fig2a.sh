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
OUTDIR=../../outputs/fig3_enjalbert2015_fig2a/
mkdir -p ${OUTDIR}
source ../../pubvenv/bin/activate

steps=10000
ibm=0.0026
minstep=0.0

psio=0.04
phio=0.04
vpsi=0.0
kpsi=0.0
vphi=0.0
kphi=0.0

### Simulate 15mM Glc condition with pFBA with and without transition
python ecDiauxie.py -b ${ibm} -n ${steps} -t 10 -x '-11.5' -s ${minstep} -p ${OUTDIR}  --run --runconsortium --ratioecgl "1.0" --runglucose -l pfba_notrans_two
python ecDiauxie.py -b ${ibm} -n ${steps} -t 10 -x '-11.5' -s ${minstep} -p ${OUTDIR}  --run --runconsortium --ratioecgl "1.0" --runglucose -l pfba_wphipsitrans_two -e '0.9' --phitransition --psitransition --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}"
### Produce plots
python ecDiauxie.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-batch_low_Glc-pfba_notrans_two_ECgl_ECac_1p0.p --ratioecgl "1.0"
python ecDiauxie.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-batch_low_Glc-pfba_wphipsitrans_two_ECgl_ECac_1p0.p --psitransition --phitransition --ratioecgl "1.0"

python agreement.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-batch_low_Glc-pfba_notrans_two_ECgl_ECac_1p0.p
python agreement.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-batch_low_Glc-pfba_wphipsitrans_two_ECgl_ECac_1p0.p


deactivate
