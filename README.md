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
First, find the path to your Python site-packages by:
```bash
python -c 'import site; print(site.getsitepackages())'
```
It is highly recommended to add the displayed path to your ```~/.bashrc```:
```bash
export PYTHON_PKGS_PATH=xxx
export GOCIA_PATH=$PYTHON_PKGS_PATH/gocia
```
Then, simply download the GOCIA tarball and untar it into the site-packages‘ path:
```bash
source ~/.bashrc
wget https://github.com/zishengz/gocia/archive/master.zip
unzip master.zip
rm -rf $GOCIA_PATH master.zip
mv gocia-master/ $GOCIA_PATH
```
After these, run the following line to test:
```bash
python -c 'import gocia'
```
If no error occurs, GOCIA should have been imported into your path!

## Tutorial

We assume the use of VASP for local optimization unless otherwise specified.

```HPC``` represents job scheduler on the cluster you use.

- SGE: Hoffman2

- SLURM: CORI

### Initial population: Structural generation

The needed files include:

- substrate.vasp
  The VASP-format structure file containing the substrate slab, with constraints.
- xxxSample.py
  Choose the structural sampling method that suit your system best.
- db2vasp.py
  Script for converting ase database files to systematically named  VASP-format files.

```bash
python xxxSample.py substrate.vasp
python db2vasp.py xxx.db
```



### Initial population: Local optimization

The needed files include:

- INCAR-1, INCAR-2, INCAR-3
  for low, mid, and high precision DFT calculations.
- KPOINTS
- init-worker.py
  The worker job that runs local optimizations.
- HPC-vasp.sh
  The shell script for submitting worker jobs.
  REMEMBER to replace the ```.bashrc``` path with yours.
- input.py
  A data file that contains the pseudo-potential path and VASP command.
- collectVASP.py
  write VASP results into a ase database file and filter out the duplicates.

```bash
chmod +x HPC-vasp.sh
for i in s0*vasp; do ./HPC-vasp.sh $i; done

python collectVASP.py
```





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

```bash
cp xxxx.db gcga.db
nohup python -u ga-HPC.py &
```



under construction 