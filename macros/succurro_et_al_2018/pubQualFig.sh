#figures/enjalbert2015_low_Glc/enjalbert2015-ecoli_core-batch_low_Glc-pfba_notrans_two_ECgl_ECac_0p95-totBM_fig2A.png
#figures/enjalbert2015_low_Glc/enjalbert2015-ecoli_core-batch_low_Glc-pfba_wphipsitrans_two_ECgl_ECac_0p95-totBM_fig2A.png
#figures/enjalbert2015_fedbatch/enjalbert2015-ecoli_core-fedbatch_low_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p95-totBM_fig6A.png
#figures/enjalbert2015_fedbatch/enjalbert2015-ecoli_core-fedbatch_high_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p75-totBM_fig6C.png
#figures/enjalbert2015_lagTime/fig4a_muExp_werr.png
#figures/enjalbert2015_lagTime/fig4b_muExp_werr.png

source ../../pubvenv/bin/activate

PQDIR=../../outputs/pqfigures
mkdir -p ${PQDIR}

### Fig. 1
#figures/enjalbert2015_oneEC/enjalbert2015-ecoli_core-batch_low_Glc-pfba_notrans_mono_ECgl-ECgl_fig2A.png
#figures/enjalbert2015_oneEC/enjalbert2015-ecoli_core-batch_low_Glc-moma_notrans_mono_ECgl-ECgl_fig2A.png

OUTDIR=../../outputs/fig1_enjalbert2015_fig2a/
python ecDiauxie.py --runsingle --runecgl  -p ${OUTDIR}  -m ${OUTDIR}/endOfSimulation-ecoli_core-batch_low_Glc-pfba_notrans_mono_ECgl.p
mv ${PQDIR}/enjalbert2015-ecoli_core-batch_low_Glc-pfba_notrans_mono_ECgl-ECgl_fig2A.png ${PQDIR}/fig1a.png
python ecDiauxie.py --runsingle --runecgl  -p ${OUTDIR}  -m ${OUTDIR}/endOfSimulation-ecoli_core-batch_low_Glc-moma_notrans_mono_ECgl.p -M
mv ${PQDIR}/enjalbert2015-ecoli_core-batch_low_Glc-moma_notrans_mono_ECgl-ECgl_fig2A.png ${PQDIR}/fig1b.png

### Fig. 2
#figures/enjalbert2015_oneEC/gene_exp_exponential_growth.png
#figures/enjalbert2015_oneEC/gene_exp_glucose_exhaustion.png

OUTDIR=../../outputs/fig2_enjalbert2015_geneexp/
rxns="ACKr PPCK PPS FBP ICL MALS PFK PYK PPC ICDHyr"
rxns0="ACKr PPCK PPS FBP ICL MALS PFK PYK"
MOMAPFBADIR=${OUTDIR}/../fig1_enjalbert2015_fig2a/

python analyseFluxes.py -i "${OUTDIR}/endOfSimulation-ecoli_core-batch_high_Ac-pfba_notrans_mono_ECgl.p ${OUTDIR}/endOfSimulation-ecoli_core-batch_low_Glc-pfba_notrans_mono_ECgl.p" -r "${rxns0}" -C -l 'Acetate Glucose' -p ${OUTDIR}
### Fig. 2b
python analyseFluxes.py -i "${MOMAPFBADIR}/endOfSimulation-ecoli_core-batch_low_Glc-moma_notrans_mono_ECgl.p ${MOMAPFBADIR}/endOfSimulation-ecoli_core-batch_low_Glc-pfba_notrans_mono_ECgl.p" -C -l 'MOMA pFBA' -p ${OUTDIR} -r "${rxns}" -G -k 0.3



### Fig 3

OUTDIR=../../outputs/fig3_enjalbert2015_fig2a/
python ecDiauxie.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-batch_low_Glc-pfba_notrans_two_ECgl_ECac_0p95.p --ratioecgl "0.95"
python ecDiauxie.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-batch_low_Glc-pfba_wphipsitrans_two_ECgl_ECac_0p95.p --psitransition --phitransition --ratioecgl "0.95"

mv ${PQDIR}/enjalbert2015-ecoli_core-batch_low_Glc-pfba_notrans_two_ECgl_ECac_0p95-totBM_fig2A.png ${PQDIR}/fig3a.png
mv ${PQDIR}/enjalbert2015-ecoli_core-batch_low_Glc-pfba_wphipsitrans_two_ECgl_ECac_0p95-totBM_fig2A.png ${PQDIR}/fig3b.png

### Fig 4
OUTDIR=../../outputs/fig4_enjalbert2015_fig6a_fig6c/

python ecDiauxie.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-fedbatch_low_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p95.p --ratioecgl "0.95" --phitransition --psitransition 
mv ${PQDIR}/enjalbert2015-ecoli_core-fedbatch_low_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p95-totBM_fig6A.png ${PQDIR}/fig4a.png

python ecDiauxie.py -p ${OUTDIR} --runconsortium -m ${OUTDIR}/endOfSimulation-ecoli_core-fedbatch_high_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p75.p --phitransition --psitransition --ratioecgl "0.75"
mv ${PQDIR}/enjalbert2015-ecoli_core-fedbatch_high_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p75-totBM_fig6C.png ${PQDIR}/fig4b.png


### Fig 5

OUTDIR=../../outputs/fig5_enjalbert2015_fig4/
psio=0.04
phio=0.04
vpsi=0.2
kpsi=30.0
vphi=0.2
kphi=5.0

python lagTimesFigure4.py -i "${OUTDIR}/endOfSimulation-ecoli_core-batch_low_Glc-pfba_wphipsitrans_two_ECgl_ECac_0p95.p" -p ${OUTDIR} -T -O --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}" -U -E

mv ${PQDIR}/fig4a_muExp_werr.png ${PQDIR}/fig5a.png

python lagTimesFigure4.py -i "${OUTDIR}/endOfSimulation-ecoli_core-batch_mixed_Ac-pfba_wphipsitrans_two_ECgl_ECac_0p75.p" -p ${OUTDIR} -T -O --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}" -U -E

mv ${PQDIR}/fig4b_muExp_werr.png ${PQDIR}/fig5b.png


deactivate
