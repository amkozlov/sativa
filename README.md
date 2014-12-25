SATIVA
======

SATIVA (**S**emi-**A**utomatic **T**axonomy **I**mprovement and **V**alidation **A**lgorithm) is a pipeline
that uses Evolutionary Placement Algorithm (EPA) to identify taxonomically mislabeled sequences
and suggest corrections 

Installation
------------

1. Make sure Python 2.6+ is installed on your system (Python 3 is not supported!)

2. Build RAxML from source by running

  $ ./install.sh


Basic usage
-----------

SATIVA requires two files as an input: alignment (FASTA or PHYLIP) and a text file with taxonomic
annotations (matched by sequence name):

  $ ./sativa.sh -s example/test.phy -t example/test.tax

Output is a text file listing the identified mislabels, confidence scores and proposed corrections.

For additional options, please refer to the online help: 

  $ ./sativa.sh -h


GUI
---

SATIVA is integrated with the most recent (unstable) version of ARB software.

Development builds: ftp://ftp.arb-silva.de/ARB/builds/

Source: http://svn.mikro.biologie.tu-muenchen.de/readonly/trunk/


Citation
--------

Manuscript is in preparation