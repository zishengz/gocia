![Picture1](./gocia_logo.png)

# G O C I A

**G**lobal **O**ptimizer for **C**lusters, **I**nterfaces, and **A**dsorbates

```GOCIA``` is a global optimization toolkit and Python modules specialized for sampling supported clusters, restructured interfaces and adsorbate configurations.

Copyright © 2020 Zisheng Zhang

[TOC]

## Requirements

- Python 3.6 or later
- ASE and its dependencies
- Natsort and LATEX (pdf report generation)

## Installation
First, install your own python environment, since HPCs usually don't give regular users write permission to the python path. To save disk space, it is recommended to install ```Miniconda```:

```shell
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh
```

Yes all the way through and ```source ~/.bashrc``` to activate the conda environment. You can run ```python``` in the terminal to check the version of the python that you are using.

Then, find the path to your Python site-packages by:

```bash
python -c 'import site; print(site.getsitepackages())'
```
It is highly recommended to add the displayed path to your ```~/.bashrc```:
```bash
export PYTHON_PKGS_PATH=xxx
export GOCIA_PATH=$PYTHON_PKGS_PATH/gocia
```
To install GOCIA, simply download the tarball and untar it into the site-packages‘ path:
```bash
source ~/.bashrc
wget https://github.com/zishengz/gocia/archive/master.zip
unzip master.zip
rm -rf $GOCIA_PATH master.zip
mv gocia-master/ $GOCIA_PATH
```
You can put the commands in the block above to a script to save time when you wanna update GOCIA the next time. After these, run the following line to test:

```bash
python -c 'import gocia'
```
If no error occurs, GOCIA should have been imported into your path!

## Tutorial

We assume the use of VASP for local optimization unless otherwise specified.

```HPC``` represents job scheduler on the cluster you use:

- slurm: CORI
- sge: Hoffman2

### 3-step local optimization

The needed files include:

- INCAR-1, INCAR-2, INCAR-3
  for low, mid, and high precision DFT calculations.
- KPOINTS
- init-worker.py
  The worker job that runs 3-step local optimizations and checks for unreasonable connectivity.
- HPC-vasp-init.sh
  The shell script for submitting worker jobs.
  REMEMBER to replace the ```.bashrc``` path with yours.
- input.py
  A data file that contains the pseudo-potential path and VASP command.

Procedure:

1. Replace the ```.bashrc``` path in ```HPC-vasp-init.sh``` with yours.
2. Put the path to your pseudo-potentials and the VASP command into ```input.py```
3. Give executable permission to the submission script by ```chmod +x HPC-vasp-init.sh```
4. Submit the job by ```./HPC-vasp-init.sh xxx.vasp```

### Initial population: Structural generation

The needed files include:

- substrate.vasp
  The VASP-format structure file containing the substrate slab, with constraints.
- xxxSample.py
  Choose the structural sampling method that suit your system best.

```bash
python xxxSample.py substrate.vasp
```

### Initial population: Local optimization

The needed files include:

- INCAR-1, INCAR-2, INCAR-3
- KPOINTS
- init-worker.py
- HPC-vasp-init.sh
- input.py
- db2vasp.py
  Script for converting ase database files to systematically named  VASP-format files.
- collectVASP.py
  write VASP results into a ase database file and filter out the duplicates.

Procudures (1-3 are the same as the 3-step opt section):

1. Replace the ```.bashrc``` path in ```HPC-vasp-init.sh``` with yours.
2. Put the path to your pseudo-potentials and the VASP command into ```input.py```
3. Give executable permission to the submission script by ```chmod +x HPC-vasp-init.sh```
4. Convert .db file to VASP-format by ```python db2vasp.py xxx.db```
5. Submit in batch by: ```for i in s0*vasp; do ./HPC-vasp-init.sh $i; done```
6. After all jobs finish, collect the results by ```python collectVASP.py```

### GCGA sampling

The needed files include:

- substrate.vasp
- INCAR-1, INCAR-2, INCAR-3
- KPOINTS
- ga-HPC.py
  The master job that runs locally (login node if permitted) and controls the job submissions.
- ga-worker.py
  The worker job that runs local optimizations and updates the population on computing nodes.
- HPC-vasp.sh
  The shell script for submitting GCGA worker jobs.
- input.py
  A data file that contains the information needed for the GCGA sampling.
- gcga.db
  The ase database file containing the initial population, obtained from the previous step.

Procedures:

1. Replace the ```.bashrc``` path in ```HPC-vasp.sh``` with yours.
2. Put the path to your pseudo-potentials, the VASP command, chemical potentials, and other GCGA parameters into ```input.py```
3. Give executable permission to the submission script by ```chmod +x HPC-vasp.sh```
4. Copy the database file from [initial population: local optimization] into ```gcga.db``` 
5. Run the GCGA master on login node by ```nohup python -u ga-HPC.py &```
6. If you want to kill the GCGA, ```touch STOP```.



Other parts are under construction...