#!/bin/bash

cpu_has_feature() {
    if [ $DARWIN -eq 1 ]; then
        SHOW="sysctl machdep.cpu.features"
    else
        SHOW="grep flags /proc/cpuinfo"
    fi
    $SHOW | grep -qi "$1"
}

USE_AVX=no
USE_AVX2=no

case `uname` in
    Darwin)
        export DARWIN=1
        CLANG_STR=`clang --version | head -1`
        if [[ $CLANG_STR =~ (LLVM ([0-9]*\.[0-9]*)svn) ]]; then
             CLANG_VERSION=${BASH_REMATCH[2]}
             COMPILER_NAME="clang $CLANG_VERSION"
             if [ "$CLANG_VERSION" \> "3.2" ]; then
                export USE_AVX=yes
             fi
             if [ "$CLANG_VERSION" \> "3.3" ]; then
                export USE_AVX2=yes
             fi
        fi
        ;;
    Linux)
        export DARWIN=0
        GCC_VERSION=`gcc -dumpversion`
        COMPILER_NAME="gcc $GCC_VERSION"
        if [ "$GCC_VERSION" \> "4.6.0" ]; then
           export USE_AVX=yes
        fi
        if [ "$GCC_VERSION" \> "4.7.0" ] || [ "$GCC_VERSION" = "4.7.0" ]; then
           export USE_AVX2=yes
        fi
        ;;
esac

if [ "$1" == "--no-avx" ]; then
  NO_AVX="yes"
  USE_AVX=no
  USE_AVX2=no
fi

if [ -z $NO_AVX ] && [ $USE_AVX != "yes" ] && cpu_has_feature avx; then
   echo "Your CPU provides AVX instuctions, but you default compiler ($COMPILER_NAME) is too old to support them."
   echo "Please consider using a more recent compiler (GCC 4.6+ or clang 3.3+) for optimal performance."
   echo "If you want to use a (slower) SSE3 RAxML version instead, please run this script with --no-avx option."
   exit
fi

if [ -z $NO_AVX ] && [ $USE_AVX2 != "yes" ] && cpu_has_feature avx2; then
   echo "Your CPU provides AVX2 instuctions, but you default compiler ($COMPILER_NAME) is too old to support them."
   echo "Please consider using a more recent compiler (GCC 4.6+ or clang 3.3+) for optimal performance."
   echo "If you want to use a (slower) SSE3 RAxML version instead, please run this script with --no-avx option."
   exit
fi

echo "Your compiler: $COMPILER_NAME"
echo "Building AVX: $USE_AVX"
echo "Building AVX2: $USE_AVX2"

cd raxml
make
cd ..
mkdir -p tmp
echo -e "\nDone!"