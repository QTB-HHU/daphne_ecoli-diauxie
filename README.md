# DAPHNE, a metabolic modeling framework to capture DynAmic Population HeterogeNEity

This modeling framework integrates Ordinary Differential Equation models with
constraint-based metabolic network models, using various Flux Balance Analysis functionalities 
implemented in COBRApy.

**The code is made available under the GNU General Public License (see LICENSE) at no warranty** and is under further development at [https://gitlab.com/asuccurro/dfba-ode-framework](https://gitlab.com/asuccurro/dfba-ode-framework). Access to the developer version can be requested by writing to *asuccurro AT gmail DOT com*.

This code corresponds to the version v1.0, used for simulating *Escherichia coli* sub-population dynamics in the work
"Emergent sub-population behavior uncovered with a community dynamic metabolic model of Escherichia coli diauxic growth" by Succurro, Segrè, Ebenhöh, 2018.

* [bioRXiv preprint](https://www.biorxiv.org/content/early/2018/10/11/291492)
* [publication]()

## Usage

Here I provide the framework code (under code/python/) and the macros (under macros/succurro_et_al_2018/)
used to obtain the Figures of the manuscript referenced above.
The code is made available under the GNU General Public License (see LICENSE) at no warranty.

Please refer to the Jupyter notebook for interactive examples to learn how to use the framework for your own analysis. A software
metapaper is in preparation.


## Run the Jupyter notebook

Load the binder environment, when asked to select the kernel choose python 2. Click on EColiSubPop.ipynb and play around.

[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/QTB-HHU/dfba-ode-framework_ecoli-diauxie/master)
	
## Installation

This code requires some external open-source packages.

We use virtual environment, if you don't have it please install like one of:

```bash
apt-get install python-virtualenv
yum install python-virtualenv
pip install virtualenv
```

### Solvers
	
Be sure to have LP solvers (default is glpk/cglpk) installed **before** pip installing cobra (you can always pip uninstall/pip install again)

In Ubuntu distributions this means e.g.:

```bash
sudo apt-get install glpk-utils
```

To install Gurobi solver (not needed) please refer to:
	
http://www.gurobi.com/documentation/6.5/quickstart_linux/software_installation_guid.html#section:Installation
	
Gurobi is installed globally by running

```bash
sudo python setup.py install
```

In virtualenv then you need to locate Gurobi path (see below)
	
### Generate a virtual environment for python2.7 like:

```bash
cd yourprojectfolder
virtualenv -p /usr/bin/python2.7 pubvenv
```

Then activate virtualenv and install the requirements:

```bash
source pubvenv/bin/activate
pip install -r req_pubvenv.txt
```

The command to leave the virtualenv is:

```bash
deactivate
```


### MOMA file modifications

The file moma.py from COBRApy has to be slightly modify to return a copy of the model:

```bash
cp filescobra/moma.py pubvenv/lib/python2.7/site-packages/cobra/flux_analysis/
```

### Gurobi installation in virtualenv

In virtualenv you need to locate Gurobi path (usually is saved in $GUROBI_HOME) and run pip like:

```bash
pip install $GUROBI_HOME
```


## Reproduce Figures from  the manuscript and the supplemental material

Please note that some simulations might be CPU intensive. More detailed information will be updated.

### Fig 1

Simulates single *E. coli* model (--runsingle --runecgl) shifting from glucose to acetate metabolism using either pFBA (default) or MOMA (-M) to solve the FBA problem.
The experimental condition is batch growth on 15 mM glucose (--runglucose).

```bash
source ./fig1_enjalbert2015_fig2a.sh
```


### Fig 2

This script uses the outputs of the previous simulation. Simulates exponential growth on glucose or acetate and obtains the flux differences as proxy for gene expression.

```bash
source ./fig2_enjalbert2015_geneexp.sh
```


### Fig 3

Simulates batch growth of two *E. coli* populations (--runconsortium) on 15 mM glucose (--runglucose), starting with 100% glucose consumers (--ratioecgl "1.0"). One simulation does not allow population shift, the other introduces noise-level population shift (--phitransition --psitransition --phioffset 0.04 --psioffset 0.04).

```bash
source ./fig3_enjalbert2015_fig2a.sh
```


### Fig 4

Simulates batch growth of two *E. coli* populations (--runconsortium) on 15 mM glucose and 4 mM acetate (--runfedlowacetate) or 15 mM glucose and 32 mM acetate (--runfedhighacetate),
starting with 100% (--ratioecgl "1.0") or 75% (--ratioecgl "0.75") glucose consumers.
Transition from one population to the other is modeled with Hill kinetics:
```bash
psio=0.04
phio=0.04
kpsi=30.0
kphi=5.0
vpsi=0.2
vphi=0.2
-e '0.9' --phitransition --psitransition --kmtransphi "${kphi}" --kmtranspsi "${kpsi}" --vmaxpsi "${vpsi}" --vmaxphi "${vphi}" --phioffset "${phio}" --psioffset "${psio}"
```

```bash
source ./fig4_enjalbert2015_fig6a_fig6c.sh
```

### Fig 5, Fig S13, Fig S16

Simulates the switch experiments. First runs mother cultures in M9G (--runglucose) and M9GA (--runmixedacetate) conditions. Then runs daughter cultures and computes lag time
using a simplified ODE model and sampling the starting population ratio from the mother cultures.

```bash
source ./fig5_enjalbert2015_fig4.sh
```

### Supplement

#### Fig S1 and S2

Simulates batch and fed-batch experiments from Varma and Palsson, 1994. Produces Rsq plots (agreement.py macro) and the flux check plots (analyseFluxes.py macro) as well.

```bash
source ./figS1-S2_varma1994_fig7_fig10.sh
```

#### Fig S3

(Accessing the simulation output files of Fig. 1) produce the flux check plots for Enjalbert et al, 2015, batch condition with glucose only.

```bash
source ./figS3_enjalbert2015.sh
``` 

#### Fig S5

```bash
source ./figS5_enjalbert_3models.sh
```

#### Fig S6

Simulate the conditions from Enjalbert et al, 2015, varying the initial biomass ratio of the two subpopulations.

```bash
source figS6_varyInitialRatio.sh
```

#### Fig S7

Simulate high acetate condition from Enjalbert et al, 2015, starting with the same initial biomass value and population ratio as the low acetate condition or as in the manuscript.

```bash
source figS7_highacetate.sh
```

#### Fig S8

```bash
source ./figS8a-d_lagtimescans.sh
source ./figS8e-h_paramscanlag_run.sh
```

#### Additional figures

```bash
source ./figS5_transitionRatesRsq.sh
```


```bash
source ./extrafig_odescans.sh
```
