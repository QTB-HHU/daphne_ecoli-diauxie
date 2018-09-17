source ../../pubvenv/bin/activate

# OPATH=../../outputs/revision_checks_figures/varma1994/
# BATCH=endOfSimulation-ecoli_core-varma_batch-varma_o2_11p5_ac_3_ECgl.p
# FEDBA=endOfSimulation-ecoli_core-varma_fedbatch-varma_o2_11p5_ac_3_ECgl.p

# python plotFluxes.py -T -p $OPATH -f ${OPATH}${BATCH} -n 'Varma, Batch growth, pFBA' -o 'varma1994_batch_ex_fluxes' -y "[['glucose_exchange', ('flux','lb')],['acetate_exchange', ('flux','lb','ub')], ['oxygen_exchange', ('flux','lb','ub')]]" -l "(7,8.5)"
# python plotFluxes.py -T -p $OPATH -f ${OPATH}${FEDBA} -n 'Varma, Fed-batch growth, pFBA' -o 'varma1994_fedbatch_ex_fluxes' -y "[['glucose_exchange', ('flux','lb')],['acetate_exchange', ('flux','lb','ub')], ['oxygen_exchange', ('flux','lb','ub')]]" -l "(7,8.5)"


OPATH=../../outputs/revision_checks_figures/fig1/
ENJPFBA=endOfSimulation-ecoli_core-batch_low_Glc-pfba_notrans_mono_ECgl.p
ENJMOMA=endOfSimulation-ecoli_core-batch_low_Glc-moma_notrans_mono_ECgl.p

python plotFluxes.py -T -p $OPATH -f ${OPATH}${ENJPFBA} -n 'Enjalbert, pFBA' -o 'enjalbert_pFBA_ex_fluxes' -y "[['glucose_exchange', ('flux','lb')],['acetate_exchange', ('flux','lb','ub')], ['oxygen_exchange', ('flux','lb','ub')]]" -l "(3.5,5)"
python plotFluxes.py -T -p $OPATH -f ${OPATH}${ENJMOMA} -n 'Enjalbert, MOMA' -o 'enjalbert_MOMA_ex_fluxes' -y "[['glucose_exchange', ('flux','lb')],['acetate_exchange', ('flux','lb','ub')], ['oxygen_exchange', ('flux','lb','ub')]]" -l "(3.5,5)"

