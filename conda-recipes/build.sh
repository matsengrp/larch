#!/bin/bash

rm -rf build
mkdir build
cd build

cmake -DCMAKE_BUILD_TYPE=${BUILD_TYPE:-"Release"} -DUSE_USHER=1 ..
make -j${NUM_THREADS:-"1"}

mkdir -p $PREFIX/lib
cp $(find . -name *.so*) $PREFIX/lib/

mkdir -p $PREFIX/bin
cp larch-usher $PREFIX/bin/larch-usher
cp dag-util $PREFIX/bin/larch-dagutil
cp dag2dot $PREFIX/bin/larch-dag2dot
cp larch-test $PREFIX/bin/larch-test
