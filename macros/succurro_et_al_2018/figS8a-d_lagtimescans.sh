#!/bin/bash 
##################################################################################
# Author: Antonella Succurro
#         a.succurro AT uni-koeln DOT de
###
# Description: reproduce Figure 5b from Succurro et al. 2018
###
# Depends on data in ../../ecoli/enjalbert2015_data/
###
# Usage: --runfedhighacetate  reproduce conditions 15mM Glc & constant ~32mM Ac  as in Enjalbert et al. 2015 Fig. 6C
#        --runfedlowacetate  reproduce conditions 15mM Glc & constant 4mM Ac  as in Enjalbert et al. 2015 Fig. 6A
#        --runconsortium --ratioecgl "0.75" runs two E. coli models with starting biomass ratio for ECgl =  0.75 TOT
#        -e '0.9' --phitransition --psitransition activates ECgl <-> ECac transitions, with efficiency 0.9
##################################################################################
OUTDIR=../../outputs/supplementary_figures/lagtime/
mkdir -p ${OUTDIR}
conda activate daphnePy2

OUTDIR=../../outputs/supplementary_figures/lagtime/muSimu/
mkdir -p ${OUTDIR}
python lagTimesFigure4.py -D -K -p ${OUTDIR}
python lagTimesFigure4.py -D -K -X -p ${OUTDIR}

OUTDIR=../../outputs/supplementary_figures/lagtime/muExp/
mkdir -p ${OUTDIR}
python lagTimesFigure4.py -D -K -U -p ${OUTDIR}
python lagTimesFigure4.py -D -K -U -M -p ${OUTDIR}
python lagTimesFigure4.py -D -K -U -N -p ${OUTDIR}

python lagTimesFigure4.py -D -K -X -U -p ${OUTDIR}
python lagTimesFigure4.py -D -K -X -U -M -p ${OUTDIR}
python lagTimesFigure4.py -D -K -X -U -N -p ${OUTDIR}


