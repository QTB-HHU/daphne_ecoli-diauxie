#!/bin/bash 
##################################################################################
# Author: Antonella Succurro
#         a.succurro AT uni-koeln DOT de
###
# Description: reproduce Figure 5 from Succurro et al. 2018
###
# Depends on data in ../../ecoli/enjalbert2015_data/
###
# Usage: --runglucose reproduce conditions 15mM Glc as in Enjalbert et al. 2015 Fig. 2A
#        --runmixedacetate reproduce conditions 15mM Glc & 32mM Ac  as in Enjalbert et al. 2015
#        --runconsortium --ratioecgl "1." runs two E. coli models with starting biomass ratio for ECgl =  0.999 TOT
#        -e '0.9' --phitransition --psitransition activates ECgl <-> ECac transitions, with efficiency 0.9
##################################################################################
OUTDIR=../../outputs/fig5_enjalbert2015_fig4/
mkdir -p ${OUTDIR}
source ../../pubvenv/bin/activate

maxsteps=10000
minstep=0.0

psio=0.04
phio=0.04
vpsi=0.2
kpsi=30.0
vphi=0.2
kphi=5.0
# psio=0.05
# phio=0.05
# vpsi=0.2
# kpsi=30.0
# vphi=0.2
# kphi=1.0

RUNMOTHER=${1:-"y"}
RUNDAUGTHER=${2:-"n"}
PLOTLAGSIMU=${3:-"n"}
PLOTLAGPAPER=${4:-"n"}
PLOTLAGODE=${5:-"y"}

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
	

case $RUNDAUGTHER in
    [Yy]* )
### Sample mother cultures grown in M9G media and simulate daughter cultures in M9G and M9A media
	OUTDIRG=${OUTDIR}/M9G/
	mkdir -p ${OUTDIRG}
	python lagTimesFigure4.py -i "${OUTDIR}/endOfSimulation-ecoli_core-batch_low_Glc-pfba_wphipsitrans_two_ECgl_ECac_0p95.p" -p ${OUTDIRG} -S --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}" > ${OUTDIRG}/runsGlc.sh
	source ${OUTDIRG}/runsGlc.sh

	### Sample mother cultures grown in M9GA media and simulate daughter cultures in M9G and M9A media
	OUTDIRGA=${OUTDIR}/M9GA/
	mkdir -p ${OUTDIRGA}
	python lagTimesFigure4.py -i "${OUTDIR}/endOfSimulation-ecoli_core-batch_mixed_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p75.p" -p ${OUTDIRGA} -S --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}" > ${OUTDIRGA}/runsAc.sh
	source ${OUTDIRGA}/runsAc.sh
	;;
esac


case $PLOTLAGSIMU in
    [Yy]* )

	### Compute the lag times and plot Fig. 6a
	x=""
	for f in ${OUTDIRG}/*p; do x="$x $f"; done
	python lagTimesFigure4.py -i "${OUTDIR}/endOfSimulation-ecoli_core-batch_low_Glc-pfba_wphipsitrans_two_ECgl_ECac_0p95.p${x}" -p ${OUTDIR} -T

	### Compute the lag times and plot Fig. 6b
	x=""
	for f in ${OUTDIRGA}/*p; do x="$x $f"; done
	python lagTimesFigure4.py -i "${OUTDIR}/endOfSimulation-ecoli_core-batch_mixed_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p75.p${x}" -p ${OUTDIR} -T
	;;
esac

case $PLOTLAGPAPER in
    [Yy]* )

	### Compute the lag times and plot Fig. 6a
	x=""
	for f in ${OUTDIRG}/*p; do x="$x $f"; done
	python lagTimesFigure4.py -i "${OUTDIR}/endOfSimulation-ecoli_core-batch_low_Glc-pfba_wphipsitrans_two_ECgl_ECac_0p95.p${x}" -p ${OUTDIR} -T -U

	### Compute the lag times and plot Fig. 6b
	x=""
	for f in ${OUTDIRGA}/*p; do x="$x $f"; done
	python lagTimesFigure4.py -i "${OUTDIR}/endOfSimulation-ecoli_core-batch_mixed_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p75.p${x}" -p ${OUTDIR} -T -U
	;;
esac

case $PLOTLAGODE in
    [Yy]* )

	python lagTimesFigure4.py -i "${OUTDIR}/endOfSimulation-ecoli_core-batch_low_Glc-pfba_wphipsitrans_two_ECgl_ECac_0p95.p" -p ${OUTDIR} -T -O --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}" -U
	python lagTimesFigure4.py -i "${OUTDIR}/endOfSimulation-ecoli_core-batch_mixed_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p75.p" -p ${OUTDIR} -T -O --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}" -U

	python lagTimesFigure4.py -i "${OUTDIR}/endOfSimulation-ecoli_core-batch_low_Glc-pfba_wphipsitrans_two_ECgl_ECac_0p95.p" -p ${OUTDIR} -T -O --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}" -U -E
	python lagTimesFigure4.py -i "${OUTDIR}/endOfSimulation-ecoli_core-batch_mixed_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p75.p" -p ${OUTDIR} -T -O --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}" -U -E

	# M9G -> stored resources
	python lagTimesFigure4.py -i "${OUTDIR}/endOfSimulation-ecoli_core-batch_low_Glc-pfba_wphipsitrans_two_ECgl_ECac_0p95.p" -p ${OUTDIR} -T -O --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}" -U -N
	python lagTimesFigure4.py -i "${OUTDIR}/endOfSimulation-ecoli_core-batch_low_Glc-pfba_wphipsitrans_two_ECgl_ECac_0p95.p" -p ${OUTDIR} -T -O --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}" -U -E -N

	# M9GA -> memory of mixed environment
	python lagTimesFigure4.py -i "${OUTDIR}/endOfSimulation-ecoli_core-batch_mixed_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p75.p" -p ${OUTDIR} -T -O --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}" -U -M
	python lagTimesFigure4.py -i "${OUTDIR}/endOfSimulation-ecoli_core-batch_mixed_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p75.p" -p ${OUTDIR} -T -O --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}" -U -E -M
	;;
esac

deactivate
