#!/bin/bash

rm -rf build
mkdir build
cd build

BUILD_TYPE=${BUILD_TYPE:-"Debug"}
USE_USHER=${USE_USHER:-"ON"}
NUM_THREADS=${NUM_THREADS:-"4"}
INCLUDE_TESTS=${INCLUDE_TESTS:-"false"}

echo "INCLUDE_TESTS: ${INCLUDE_TESTS}"
echo "BUILD_TYPE: ${BUILD_TYPE}"
echo "USE_USHER: ${USE_USHER}"
echo "NUM_THREADS: ${NUM_THREADS}"

cmake -DCMAKE_BUILD_TYPE=${BUILD_TYPE} -DUSE_USHER=${USE_USHER} ..
make -j${NUM_THREADS}

mkdir -p $PREFIX/lib
cp $(find . -name *.so*) $PREFIX/lib/

mkdir -p $PREFIX/bin
cp larch-usher $PREFIX/bin/larch-usher
cp larch-dagutil $PREFIX/bin/larch-dagutil
cp larch-dag2dot $PREFIX/bin/larch-dag2dot

# if [[ ${INCLUDE_TESTS} == true ]]; then
    cp larch-test $PREFIX/bin/larch-test
# fi

