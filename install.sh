#!/bin/bash

cpu_has_feature() {
    if [ $DARWIN -eq 1 ]; then
        SHOW="sysctl machdep.cpu.features"
    else
        SHOW="grep flags /proc/cpuinfo"
    fi
    $SHOW | grep -qi "$1"
}

case `uname` in
    Darwin)
        export DARWIN=1
        ;;
    Linux)
        export DARWIN=0
        ;;
esac

if cpu_has_feature avx; then
    export USE_AVX=yes
fi

cd raxml
make
cd ..
echo "\nDone!"