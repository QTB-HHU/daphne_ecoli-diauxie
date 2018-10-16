#!/bin/bash 
##################################################################################
# Author: Antonella Succurro
#         a.succurro AT uni-koeln DOT de
###
# Description: Figure 3 and 4 testing different initial ratios from Succurro et al. 2018
###
# Depends on data in ../../ecoli/enjalbert2015_data/
###
# Usage: --runglucose reproduce conditions 15mM Glc as in Enjalbert et al. 2015 Fig. 2A
#        --runfedlowacetate  reproduce conditions 15mM Glc & constant 4mM Ac  as in Enjalbert et al. 2015 Fig. 6A
#        --runconsortium --ratioecgl "1." runs two E. coli models with starting biomass ratio for ECgl =  0.999 TOT
#        -e '0.9' --phitransition --psitransition activates ECgl <-> ECac transitions, with efficiency 0.9
##################################################################################
OUTDIR=../../outputs/vary_initial_ratio/
mkdir -p ${OUTDIR}
source ../../pubvenv/bin/activate

c1="95 75 5 25 05"
array1=( $c1 )
c2="0p95 0p75 0p5 0p25 0p05"
array2=( $c2 )

for ((i=0;i<5;i++)); do
    rv="0.${array1[$i]}"
    rs="${array2[$i]}"
    echo "${rv} ${rs}"

    ########### Figure 3
    steps=10000
    ibm=0.0027
    minstep=0.0

    psio=0.04
    phio=0.04
    vpsi=0.0
    kpsi=0.0
    vphi=0.0
    kphi=0.0	

    ### Simulate 15mM Glc condition with pFBA with and without transition
    python ecDiauxie.py -b ${ibm} -n ${steps} -t 10 -x '-11.5' -s ${minstep} -p ${OUTDIR}  --run --runconsortium --ratioecgl "${rv}" --runglucose -l pfba_notrans_two
    python ecDiauxie.py -b ${ibm} -n ${steps} -t 10 -x '-11.5' -s ${minstep} -p ${OUTDIR}  --run --runconsortium --ratioecgl "${rv}" --runglucose -l pfba_wphipsitrans_two -e '0.9' --phitransition --psitransition --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}"
    ### Produce plots
    #python ecDiauxie.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-batch_low_Glc-pfba_notrans_two_ECgl_ECac_${rs}.p --ratioecgl "${rv}"
    #python ecDiauxie.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-batch_low_Glc-pfba_wphipsitrans_two_ECgl_ECac_${rs}.p --psitransition --phitransition --ratioecgl "${rv}"

    #python agreement.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-batch_low_Glc-pfba_notrans_two_ECgl_ECac_${rs}.p
    #python agreement.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-batch_low_Glc-pfba_wphipsitrans_two_ECgl_ECac_${rs}.p

    # Figure 4
    maxsteps=10000
    ibm=0.0038
    minstep=0.0
    
    psio=0.04
    phio=0.04
    kpsi=30.0
    kphi=5.0
    vpsi=0.2
    vphi=0.2

    python ecDiauxie.py -b ${ibm} -n ${maxsteps} -t 13 -x '-11.5' -s ${minstep} -p ${OUTDIR} --run --runconsortium --ratioecgl "${rv}" --runfedlowacetate -l pfba_wphipsitrans_two -e '0.9' --phitransition --psitransition --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}"
    #python ecDiauxie.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-fedbatch_low_Ac-pfba_wphipsitrans_two_ECgl_ECac_${rs}.p --ratioecgl "${rv}" --phitransition --psitransition 

    ibm=0.006

    python ecDiauxie.py -b ${ibm} -n ${maxsteps} -t 11 -x '-11.5' -s ${minstep} -p ${OUTDIR} --run --runconsortium --ratioecgl "${rv}" --runfedhighacetate -l pfba_wphipsitrans_two -e '0.9' --phitransition --psitransition --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}"
    #python ecDiauxie.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-fedbatch_high_Ac-pfba_wphipsitrans_two_ECgl_ECac_${rs}.p --phitransition --psitransition --ratioecgl "${rv}"

done

python plotVaryRatio.py -i '../../outputs/vary_initial_ratio/endOfSimulation-ecoli_core-fedbatch_high_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p05.p ../../outputs/vary_initial_ratio/endOfSimulation-ecoli_core-fedbatch_high_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p25.p ../../outputs/vary_initial_ratio/endOfSimulation-ecoli_core-fedbatch_high_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p5.p ../../outputs/vary_initial_ratio/endOfSimulation-ecoli_core-fedbatch_high_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p75.p ../../outputs/vary_initial_ratio/endOfSimulation-ecoli_core-fedbatch_high_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p95.p' -p '../../outputs/vary_initial_ratio/'

python plotVaryRatio.py -i '../../outputs/vary_initial_ratio/endOfSimulation-ecoli_core-batch_low_Glc-pfba_wphipsitrans_two_ECgl_ECac_0p05.p ../../outputs/vary_initial_ratio/endOfSimulation-ecoli_core-batch_low_Glc-pfba_wphipsitrans_two_ECgl_ECac_0p25.p ../../outputs/vary_initial_ratio/endOfSimulation-ecoli_core-batch_low_Glc-pfba_wphipsitrans_two_ECgl_ECac_0p5.p ../../outputs/vary_initial_ratio/endOfSimulation-ecoli_core-batch_low_Glc-pfba_wphipsitrans_two_ECgl_ECac_0p75.p ../../outputs/vary_initial_ratio/endOfSimulation-ecoli_core-batch_low_Glc-pfba_wphipsitrans_two_ECgl_ECac_0p95.p' -p '../../outputs/vary_initial_ratio/'

python plotVaryRatio.py -i '../../outputs/vary_initial_ratio/endOfSimulation-ecoli_core-fedbatch_low_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p05.p ../../outputs/vary_initial_ratio/endOfSimulation-ecoli_core-fedbatch_low_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p25.p ../../outputs/vary_initial_ratio/endOfSimulation-ecoli_core-fedbatch_low_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p5.p ../../outputs/vary_initial_ratio/endOfSimulation-ecoli_core-fedbatch_low_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p75.p ../../outputs/vary_initial_ratio/endOfSimulation-ecoli_core-fedbatch_low_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p95.p' -p '../../outputs/vary_initial_ratio/'


deactivate
