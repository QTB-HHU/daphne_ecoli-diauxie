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


### Fig 1

Simulates single *E. coli* model (--runsingle --runecgl) shifting from glucose to acetate metabolism using either pFBA (default) or MOMA (-M) to solve the FBA problem.
The experimental condition is batch growth on 15 mM glucose (--runglucose).

```bash
source ./fig1_enjalbert2015_fig2a.sh
```

#### time output

```bash
/usr/bin/time -v ./fig1_enjalbert2015_fig2a.sh

User time (seconds): 964.89
System time (seconds): 15.70
Percent of CPU this job got: 98%
Elapsed (wall clock) time (h:mm:ss or m:ss): 16:36.93
Average shared text size (kbytes): 0
Average unshared data size (kbytes): 0
Average stack size (kbytes): 0
Average total size (kbytes): 0
Maximum resident set size (kbytes): 1603088
Average resident set size (kbytes): 0
Major (requiring I/O) page faults: 0
Minor (reclaiming a frame) page faults: 8262327
Voluntary context switches: 12100
Involuntary context switches: 2094534
Swaps: 0
File system inputs: 0
File system outputs: 484304
Socket messages sent: 0
Socket messages received: 0
Signals delivered: 0
Page size (bytes): 4096
Exit status: 0
```

### Fig 2

This script uses the outputs of the previous simulation. Simulates exponential growth on glucose or acetate and obtains the flux differences as proxy for gene expression.

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
Average shared text size (kbytes): 0
Average unshared data size (kbytes): 0
Average stack size (kbytes): 0
Average total size (kbytes): 0
Maximum resident set size (kbytes): 251816
Average resident set size (kbytes): 0
Major (requiring I/O) page faults: 0
Minor (reclaiming a frame) page faults: 455864
Voluntary context switches: 1786
Involuntary context switches: 1292846
Swaps: 0
File system inputs: 0
File system outputs: 2472
Socket messages sent: 0
Socket messages received: 0
Signals delivered: 0
Page size (bytes): 4096
Exit status: 0
```

### Fig 3

Simulates batch growth of two *E. coli* populations (--runconsortium) on 15 mM glucose (--runglucose), starting with 100% glucose consumers (--ratioecgl "1.0"). One simulation does not allow population shift, the other introduces noise-level population shift (--phitransition --psitransition --phioffset 0.04 --psioffset 0.04).

```bash
source ./fig3_enjalbert2015_fig2a.sh
```

#### time output

```bash
/usr/bin/time -v ./fig3_enjalbert2015_fig2a.sh

User time (seconds): 76.10
System time (seconds): 9.32
Percent of CPU this job got: 111%
Elapsed (wall clock) time (h:mm:ss or m:ss): 1:16.36
Average shared text size (kbytes): 0
Average unshared data size (kbytes): 0
Average stack size (kbytes): 0
Average total size (kbytes): 0
Maximum resident set size (kbytes): 239676
Average resident set size (kbytes): 0
Major (requiring I/O) page faults: 0
Minor (reclaiming a frame) page faults: 809509
Voluntary context switches: 2397
Involuntary context switches: 2146733
Swaps: 0
File system inputs: 96
File system outputs: 63432
Socket messages sent: 0
Socket messages received: 0
Signals delivered: 0
Page size (bytes): 4096
Exit status: 0
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
#### time output

```bash
/usr/bin/time -v ./fig4_enjalbert2015_fig6a_fig6c.sh

Command being timed: "./fig4_enjalbert2015_fig6a_fig6c.sh"
User time (seconds): 65.88
System time (seconds): 5.64
Percent of CPU this job got: 107%
Elapsed (wall clock) time (h:mm:ss or m:ss): 1:06.71
Average shared text size (kbytes): 0
Average unshared data size (kbytes): 0
Average stack size (kbytes): 0
Average total size (kbytes): 0
Maximum resident set size (kbytes): 214688
Average resident set size (kbytes): 0
Major (requiring I/O) page faults: 0
Minor (reclaiming a frame) page faults: 612717
Voluntary context switches: 1283
Involuntary context switches: 1021841
Swaps: 0
File system inputs: 0
File system outputs: 59928
Socket messages sent: 0
Socket messages received: 0
Signals delivered: 0
Page size (bytes): 4096
Exit status: 0
```

### Fig 5

Simulates the switch experiments. First runs mother cultures in M9G (--runglucose) and M9GA (--runmixedacetate) conditions. Then runs daughter cultures and computes lag time
using a simplified ODE model and sampling the starting population ratio from the mother cultures.

```bash
source ./fig5_enjalbert2015_fig4.sh
```

#### time output

```bash
/usr/bin/time -v ./fig5_enjalbert2015_fig4.sh

User time (seconds): 165.74
System time (seconds): 66.93
Percent of CPU this job got: 150%
Elapsed (wall clock) time (h:mm:ss or m:ss): 2:35.03
Average shared text size (kbytes): 0
Average unshared data size (kbytes): 0
Average stack size (kbytes): 0
Average total size (kbytes): 0
Maximum resident set size (kbytes): 616420
Average resident set size (kbytes): 0
Major (requiring I/O) page faults: 0
Minor (reclaiming a frame) page faults: 2520492
Voluntary context switches: 14582
Involuntary context switches: 4792230
Swaps: 0
File system inputs: 0
File system outputs: 80528
Socket messages sent: 0
Socket messages received: 0
Signals delivered: 0
Page size (bytes): 4096
Exit status: 0
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
