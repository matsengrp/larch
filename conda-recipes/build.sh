#!/bin/bash

echo "Current directory: $(pwd)"
ls -l

rm -rf build

mkdir -p build
cd build
echo "Build directory contents:"
ls -l

cmake -DCMAKE_BUILD_TYPE=Release -DUSE_USHER=1 ..
make

mkdir -p $PREFIX/lib
cp $(find . -name *.so*) $PREFIX/lib/

mkdir -p $PREFIX/bin
cp larch-usher $PREFIX/bin/larch-usher
cp dag-util $PREFIX/bin/larch-dagutil
cp dag2dot $PREFIX/bin/larch-dag2dot
cp larch-test $PREFIX/bin/larch-test
