for vpsi in 0.2 0.5 0.8; 
do
    for vphi in 0.2 0.5 0.8;
    do
	for kpsi in 10.0 30.0;
	do
	    for kphi in 5.0;
	    do
		echo "${vpsi} ${vphi} ${kpsi} ${kphi}"
		source figs15_paramscanlag_run.sh y y "${vpsi}" "${vphi}" "${kpsi}" "${kphi}"
	    done
	done
    done
done

OUTDIR=../../outputs/supplementary_figures/lagtime/scans/
mkdir -p ${OUTDIR}
source ../../pubvenv/bin/activate

mkdir -p ${OUTDIR}/kpsi10/
mkdir -p ${OUTDIR}/kpsi30/

python ratioFigureSup.py -m 'M9G' -p ${OUTDIR} --kmtranspsi 10.0
python ratioFigureSup.py -m 'M9G' -p ${OUTDIR} --kmtranspsi 30.0
python ratioFigureSup.py -m 'M9GA' -p ${OUTDIR} --kmtranspsi 10.0
python ratioFigureSup.py -m 'M9GA' -p ${OUTDIR} --kmtranspsi 30.0
