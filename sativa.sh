#!/bin/sh

# Usage info
show_help() {
cat << EOF
Usage: ${0##*/} -t TAXFILE -s SEQFILE [-n NAME] [-o OUTDIR] [-T THREADS] [-d LEVEL] [-h]
Find taxonomically mislabeled sequences

    -t TAXFILE  Taxonomy file.
    -s SEQFILE  Alignment file (FASTA or PHYLIP). Sequence IDs must match those in TAXFILE.
    -n NAME     run name. Default: identical to TAXFILE
    -T THREADS  number of threads. Default: equals to number of CPU cores
    -o OUTDIR   directory for output files. Default: same directory as TAXFILE
    -d LEVEL    debug level: 0, 1 or 2. Default: 0 (debug off).
    -h          display this help and exit
EOF
}

cpu_get_cores() {
    case `uname` in
        Darwin)
            sysctl -n hw.ncpu
            ;;
        Linux)
            grep -c "^processor" /proc/cpuinfo
            ;;
    esac
}

# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

MODE=thorough

EPACDIR=epac
TMPDIR=tmp

CFG=sativa.cfg

DEBUG=0

while getopts "hd:t:s:n:o:T:" opt; do
    case "$opt" in
    h)
        show_help
        exit 0
        ;;
    d)  DEBUG=$OPTARG
        ;;
    t)  INFILE_TAX=$OPTARG
        ;;
    s)  INFILE_ALI=$OPTARG
        ;;
    n)  NAME=$OPTARG
        ;;
    o)  OUTDIR=$OPTARG
        ;;
    T)  THREADS=$OPTARG
        ;;
    esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

if [ -z "$INFILE_TAX" ] || [ -z "$INFILE_ALI" ]; then
  show_help
  exit 0
fi

if [ -z "$OUTDIR" ]; then
  OUTDIR=$(dirname "$INFILE_TAX")
fi

if [ -z "$NAME" ]; then
  NAME=$(basename "$INFILE_TAX")
fi

if [ -z "$THREADS" ]; then
  THREADS=`cpu_get_cores`
fi

JSON=$TMPDIR/$NAME.json

if [ "$DEBUG" -eq 2 ]; then
  DBG="-debug"
else
  DBG=""
fi

$EPACDIR/epa_trainer.py -t $INFILE_TAX -s $INFILE_ALI -r $JSON -c $CFG -v -C -no-hmmer -tmpdir $TMPDIR -m $MODE -T $THREADS $DBG

$EPACDIR/find_mislabels.py -r $JSON -c $CFG -o $OUTDIR -n $NAME -tmpdir $TMPDIR -T $THREADS $DBG

L1OUT_JPLACE=$TMPDIR/RAxML_leaveOneOutResults.l1out_seq_${NAME}.jplace
FINAL_JPLACE=$TMPDIR/RAxML_portableTree.final_epa_${NAME}.jplace
if [ "$DEBUG" -ge 1 ]; then
  mv $JSON $OUTDIR

  mv $L1OUT_JPLACE $OUTDIR/${NAME}_l1out.jplace

  mv $FINAL_JPLACE $OUTDIR/${NAME}_final.jplace
else
  rm $JSON $L1OUT_JPLACE $FINAL_JPLACE
fi
