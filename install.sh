#!/bin/bash

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



if cpu_has_feature avx; then
    export USE_AVX=yes
fi

cd raxml
make
cd ..
echo "Done!"