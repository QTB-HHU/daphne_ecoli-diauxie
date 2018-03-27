# Modeling framework for dynamic metabolic modeling of ecosystems

This modeling framework integrates Ordinary Differential Equation models with
constraint-based metabolic network models, using various Flux Balance Analysis functionalities 
implemented in COBRApy.

**The code is made available under the GNU General Public License (see LICENSE) at no warranty** and is under further development at [https://gitlab.com/asuccurro/dfba-ode-framework](https://gitlab.com/asuccurro/dfba-ode-framework). Access to the developer version can be requested by writing to *asuccurro AT gmail DOT com*.

This code corresponds to the version v1.0, used for simulating *Escherichia coli* sub-population dynamics in the work
"Emerging sub-population dynamics uncovered with a community model of Escherichia coli diauxic growth" by Succurro, Segrè, Ebenhöh, 2018.

* [bioRXiv preprint]()
* [publication]()


## Usage

We provide the framework code (under code/python/) and the macros (under macros/succurro_et_al_2018/)
used to obtain the Figures of the manuscript and the supplemental material.
The code is made available under the GNU General Public License (see LICENSE) at no warranty.
Please note that some of the simulations are rather computationally intensive.

### Fig 1

```bash
source ./fig1_enjalbert2015_fig2a.sh
```

#### time output

```bash
/usr/bin/time -v ./fig1_enjalbert2015_fig2a.sh
```

### Fig 2

```bash
source ./fig2_enjalbert2015_geneexp.sh
```

#### time output

```bash
/usr/bin/time -v ./fig2_enjalbert2015_geneexp.sh

Command being timed: "./fig2_enjalbert2015_geneexp.sh"
User time (seconds): 17.94
System time (seconds): 8.84
Percent of CPU this job got: 155%
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:17.21
```

### Fig 3

```bash
source ./fig3_enjalbert2015_fig2a.sh
```

#### time output

```bash
/usr/bin/time -v ./fig3_enjalbert2015_fig2a.sh
```



### Fig 4

```bash
source ./fig4_enjalbert2015_fig6a_fig6c.sh
```
#### time output

```bash
/usr/bin/time -v ./fig4_enjalbert2015_fig6a_fig6c.sh

Command being timed: "./fig4_enjalbert2015_fig6a_fig6c.sh"
User time (seconds): 65.88
System time (seconds): 5.64
Percent of CPU this job got: 107%
Elapsed (wall clock) time (h:mm:ss or m:ss): 1:06.71
```

### Fig 5

```bash
source ./fig5_enjalbert2015_fig4.sh
```

#### time output

```bash
/usr/bin/time -v
 ./fig5_enjalbert2015_fig4.sh
```


### Supplement

```bash
source ./figs12_lagtimescans.sh
```

```bash
source ./figs13_paramscanlag.sh
```

```bash
source ./figs1_varma1994_fig7_fig10.sh
```

```bash
source ./figs5-s7_odescans.sh
```

```bash
source ./figs9-s11_transitionrates.sh
```


				 
				 



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
