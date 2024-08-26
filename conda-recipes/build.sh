#!/bin/bash

mkdir build
cd build
cmake -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_BUILD_TYPE=Debug -DUSE_USHER=1 ..
make -j${CPU_COUNT}

mkdir -p $PREFIX/lib
cp $(find . -name *.so*) $PREFIX/lib/

mkdir -p $PREFIX/bin
cp larch-usher $PREFIX/bin/
cp larch-dagutil $PREFIX/bin/
cp larch-dag2dot $PREFIX/bin/
cp larch-test $PREFIX/bin/
