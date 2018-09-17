OUTDIR=../../outputs/revision_checks_figures/fig1/
mkdir -p ${OUTDIR}
source ../../pubvenv/bin/activate

maxsteps=10000
ibm=0.0027

### Simulate 15mM Glc condition with MOMA and pFBA
for a in 3 4 5 6 7 8 9 10;
    do
	sta=${a//\./p}
	python ecDiauxie_revision.py -b ${ibm} -n ${maxsteps} -t 10. -x '-11.5' -q "${a}" -p ${OUTDIR} --run --runsingle --runecgl --ratioecgl "1." --runglucose -l pfba_ac_upt_${sta}

	python ecDiauxie_revision.py --runsingle --runecgl  -p ${OUTDIR}  -m ${OUTDIR}/endOfSimulation-ecoli_core-batch_low_Glc-pfba_ac_upt_${sta}_ECgl.p
done
