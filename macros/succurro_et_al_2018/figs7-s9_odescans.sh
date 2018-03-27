#!/bin/bash 
##################################################################################
# Author: Antonella Succurro
#         a.succurro AT uni-koeln DOT de
###
# Description: 
###
# Depends on data in ../../ecoli/enjalbert2015_data/
###
# Usage: --runglucose reproduce conditions 15mM Glc as in Enjalbert et al. 2015 Fig. 2A
#        --runfedlowacetate  reproduce conditions 15mM Glc & constant 4mM Ac  as in Enjalbert et al. 2015 Fig. 6A
#        --runsingle --runecgl --ratioecgl "1." runs a single E. coli model
#        -M runs MOMA, else pFBA
##################################################################################
OUTDIR=../../outputs/supplementary_figures/
mkdir -p ${OUTDIR}
source ../../pubvenv/bin/activate

maxsteps=10000
ibm=0.08

### Supp Mat ODE sys
ODIR=${OUTDIR}/gamma/
mkdir -p $ODIR

python odeSolRatioBM.py -p $ODIR -O -A -r 1.0
python odeSolRatioBM.py -p $ODIR -O -G -r 1.0
python odeSolRatioBM.py -p $ODIR -O -D -r 1.0
python odeSolRatioBM.py -p $ODIR -O -A -r 0.5
python odeSolRatioBM.py -p $ODIR -O -G -r 0.5
python odeSolRatioBM.py -p $ODIR -O -D -r 0.5
python odeSolRatioBM.py -p $ODIR -O -A -r 0.0
python odeSolRatioBM.py -p $ODIR -O -G -r 0.0
python odeSolRatioBM.py -p $ODIR -O -D -r 0.0

python odeSolRatioBM.py -p $ODIR -K

deactivate
