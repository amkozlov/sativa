SATIVA
======

SATIVA (**S**emi-**A**utomatic **T**axonomy **I**mprovement and **V**alidation **A**lgorithm) is a pipeline
that uses Evolutionary Placement Algorithm (EPA) to identify taxonomically mislabeled sequences
and suggest corrections 

Installation
------------

Currently, only Linux and OSX (Mac) systems are supported. 

1. Make sure Python 2.6+ is installed (Python 3 is not supported!)

2. Make sure you have a recent C compiler (we recommend GCC 4.6+ / clang 3.3+ for AVX support).
   If you have an up-to-date OS distribution (Ubuntu 12.04+, OSX 10.8+ etc.), there is nothing to worry about.
   In a cluster environment, you might need to select an appropriate compiler version, e.g.:

   `module load gcc/4.7.0`

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

Sample command line:

```cd example
   ../sativa.py -s test.phy -t test.tax -x BAC
```

Output is a text file listing the identified mislabels, confidence scores and proposed corrections.

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

Alexey M. Kozlov, Jiajie Zhang, Pelin Yilmaz, Frank Oliver Glöckner and Alexandros Stamatakis.
**Phylogeny-aware Identification and Correction of Taxonomically Mislabeled Sequences.** *In preparation.*