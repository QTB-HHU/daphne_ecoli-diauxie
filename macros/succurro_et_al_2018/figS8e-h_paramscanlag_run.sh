#!/bin/bash 
##################################################################################
# Author: Antonella Succurro
#         a.succurro AT uni-koeln DOT de
###
# Description: 
###
# Depends on data in ../../ecoli/enjalbert2015_data/
###
# Usage: 
##################################################################################
OUTDIR=../../outputs/supplementary_figures/lagtime/scans/
mkdir -p ${OUTDIR}
source ../../pubvenv/bin/activate

maxsteps=10000
minstep=0.0

psio=0.04
phio=0.04

######################
RUNMOTHER=${1:-"n"}
PLOTLAGODE=${2:-"n"}
vpsi=${3:-"0.2"}
vphi=${4:-"0.2"}
kpsi=${5:-"30.0"}
kphi=${6:-"1.0"}

OUTDIR=../../outputs/supplementary_figures/lagtime/scans/vpsi${vpsi}_kpsi${kpsi}_vphi${vphi}_kphi${kphi}/
mkdir -p ${OUTDIR}

case $RUNMOTHER in
    [Yy]* )
	### Simulate mother cultures
	# Glucose batch, like Enjalbert Fig2A
	ibm=0.0027
	python ecDiauxie.py -b ${ibm} -n ${maxsteps} -t 10 -x '-11.3' -s ${minstep} -p ${OUTDIR} --run --runconsortium --ratioecgl "0.95" --runglucose -l pfba_wphipsitrans_two -e '0.9' --phitransition --psitransition --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}"
	python ecDiauxie.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-batch_low_Glc-pfba_wphipsitrans_two_ECgl_ECac_0p95.p --psitransition --phitransition --ratioecgl "0.95"
	# Glucose+Acetate, similar to Enjalbert Fig6C
	ibm=0.006
	python ecDiauxie.py -b ${ibm} -n ${maxsteps} -t 10 -x '-11.3' -s ${minstep} -p ${OUTDIR} --run --runconsortium --ratioecgl "0.75" --runmixedacetate -l pfba_wphipsitrans_two -e '0.9' --phitransition --psitransition --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}"
	python ecDiauxie.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-batch_mixed_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p75.p --psitransition --phitransition --ratioecgl "0.75"

	;;
esac

#px=0p95
px=0p75

case $PLOTLAGODE in
    [Yy]* )

	python lagTimesFigure4.py -i "${OUTDIR}/endOfSimulation-ecoli_core-batch_low_Glc-pfba_wphipsitrans_two_ECgl_ECac_0p95.p" -p ${OUTDIR} -T -O --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}" -U
	python lagTimesFigure4.py -i "${OUTDIR}/endOfSimulation-ecoli_core-batch_mixed_Ac-pfba_wphipsitrans_two_ECgl_ECac_${px}.p" -p ${OUTDIR} -T -O --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}" -U

	;;
    
    [Ss]* )

	python lagTimesFigure4.py -i "${OUTDIR}/endOfSimulation-ecoli_core-batch_low_Glc-pfba_wphipsitrans_two_ECgl_ECac_0p95.p" -p ${OUTDIR} -T -O --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}"
	python lagTimesFigure4.py -i "${OUTDIR}/endOfSimulation-ecoli_core-batch_mixed_Ac-pfba_wphipsitrans_two_ECgl_ECac_${px}.p" -p ${OUTDIR} -T -O --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}"

	;;
esac

deactivate
