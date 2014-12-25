#!/bin/bash
set -e

# return true if argument is file in path executable by user
can_run() {
    which "$1" && test -x `which "$1"`
}

# show error in ARB and exit
report_error() {
    SELF=`basename "$0"`
    echo "$SELF failed with: $*"
}

cpu_has_feature() {
    case `uname` in
        Darwin)
            SHOW="sysctl machdep.cpu.features"
            ;;
        Linux)
            SHOW="grep flags /proc/cpuinfo"
            ;;
    esac
    $SHOW | grep -qi "$1"
}

###### main #####

# select RAxML binary
BASEDIR=$(dirname $0)
RAXML_STEM="${BASEDIR}/raxmlHPC8"
if cpu_has_feature avx && can_run ${RAXML_STEM}-AVX.PTHREADS; then
    RAXML=${RAXML_STEM}-AVX.PTHREADS
elif cpu_has_feature sse3 && can_run ${RAXML_STEM}-SSE3.PTHREADS; then
    RAXML=${RAXML_STEM}-SSE3.PTHREADS
elif can_run ${RAXML_STEM}-PTHREADS; then
    RAXML=${RAXML_STEM}-PTHREADS
else
    report_error Could not find working RAxML binary. 
    exit -1
fi

$RAXML $@
