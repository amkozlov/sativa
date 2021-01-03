SATIVA
======

[![Build Status](https://travis-ci.org/amkozlov/sativa.svg?branch=master)](https://travis-ci.org/amkozlov/sativa)
[![License](https://img.shields.io/badge/license-GPL3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0.en.html)

SATIVA (**S**emi-**A**utomatic **T**axonomy **I**mprovement and **V**alidation **A**lgorithm) is a pipeline
that uses Evolutionary Placement Algorithm (EPA, [1]) to identify taxonomically mislabeled sequences
and suggest corrections. Internally, SATIVA relies on RAxML [2] for likelihood computations as well as on 
the ETE library[3] for tree topology manupulations in Python.

Installation
------------

Currently, only Linux and OSX (Mac) systems are supported. 

1. Make sure Python 3 is installed (Python 2 is not supported!)

2. Make sure you have a recent C compiler (we recommend GCC 4.6+ / clang 3.3+ for AVX support).
   If you have an up-to-date OS distribution (Ubuntu 12.04+, OSX 10.8+ etc.), there is nothing to worry about.
   In a cluster environment, you might need to select an appropriate compiler version, e.g.:

   `module load gcc`

   (please refer to your cluster documentation for details)

3. Run the installation script

  `./install.sh`

  If you are getting compilation errors, try to disable AVX:

  `./install.sh --no-avx`

Basic usage
-----------

SATIVA requires two files as an input: alignment (FASTA or PHYLIP) and a text file with taxonomic
annotations (matched by sequence name). Furtermore, you must choose the nomenclature code via the
-x option (e.g., BAC(teriological) for Bacteria and Archaea).

Sample command line to run SATIVA with 2 threads:

```sh
   cd example
   ../sativa.py -s test.phy -t test.tax -x BAC -T 2
```

Output is a text file which contains a list of identified mislabels, along with the corresponding 
confidence scores and proposed taxonomic corrections.

**Parallelization note**: If you omit the `-T` parameter, SATIVA will start one thread per each logical CPU
in your system. Although this is usually what you want, it might lead to a *major* slowdown
if some of the CPUs are already reserved by other running programs (e.g., if you run SATIVA on 
a shared server). If you encounter this problem, please try reducing the number of threads with `-T`!

**Handling non-preferred synonyms**: You can use `-Y` parameter to specify a file with the list of
equivalent name groups (synonyms), e.g.:

```sh
   cd example
   ../sativa.py -s test.phy -t test.tax -x BAC -T 2 -n syntest -Y synonym.txt
```

First name of each group will be considered the preferred synonym (primary name), and 
will be used in place of all other (synonymous) names in the group. 
An example synonym definition can be found in `synonym.txt` file.

For additional options, please refer to the online help:

  `./sativa.py -h`

GUI
---

SATIVA is integrated with the most recent (unstable) version of ARB software.

Development builds: ftp://ftp.arb-silva.de/ARB/builds/

Source: http://svn.mikro.biologie.tu-muenchen.de/readonly/trunk/


Support
-------

For the time being, please direct your questions to the RAxML google group:

https://groups.google.com/forum/?hl=en#!forum/raxml


Citation
--------

Alexey M. Kozlov, Jiajie Zhang, Pelin Yilmaz, Frank Oliver Glöckner and Alexandros Stamatakis (2016)
**Phylogeny-aware Identification and Correction of Taxonomically Mislabeled Sequences.**
*Nucleic Acids Research* [open access](https://nar.oxfordjournals.org/content/early/2016/05/09/nar.gkw396.full)


References
----------

[1] Berger, S. A., Krompass, D., and Stamatakis, A. (2011) 
**Performance, Accuracy, and Web Server for Evolutionary Placement of Short Sequence Reads under Maximum Likelihood.**
*Systematic Biology*, 60(3), 291–302.
doi:[10.1093/sysbio/syr010](http://sysbio.oxfordjournals.org/content/60/3/291)

[2] Stamatakis A. (2014)
**RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies.**
*Bioinformatics*, 30(9): 1312-1313.
doi:[10.1093/bioinformatics/btu033](http://dx.doi.org/10.1093/bioinformatics/btu033)

[3] Huerta-Cepas, J., Dopazo, J., and Gabaldon, T. (2010) 
**ETE: a python Environment for Tree Exploration.**
*BMC bioinformatics*, 11(1), 24.
doi:[10.1186/1471-2105-11-24](http://www.biomedcentral.com/1471-2105/11/24)
