#!/bin/bash 
OUTDIR=../../outputs/supplementary_figures/high_acetate/
mkdir -p ${OUTDIR}
conda activate daphnePy2

maxsteps=10000
ibm=0.0038
minstep=0.0

psio=0.04
phio=0.04
kpsi=30.0
kphi=5.0
vpsi=0.2
vphi=0.2

python ecDiauxie.py -b ${ibm} -n ${maxsteps} -t 11 -x '-11.5' -s ${minstep} -p ${OUTDIR} --run --runconsortium --ratioecgl "0.95" --runfedhighacetate -l pfba_wphipsitrans_two -e '0.9' --phitransition --psitransition --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}"
python ecDiauxie.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-fedbatch_high_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p95.p --phitransition --psitransition --ratioecgl "0.95" --shiftgetime

ibm=0.006

python ecDiauxie.py -b ${ibm} -n ${maxsteps} -t 11 -x '-11.5' -s ${minstep} -p ${OUTDIR} --run --runconsortium --ratioecgl "0.75" --runfedhighacetate -l pfba_wphipsitrans_two -e '0.9' --phitransition --psitransition --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}"
python ecDiauxie.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-fedbatch_high_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p75.p --phitransition --psitransition --ratioecgl "0.75"


