LARCHBUILDPATH="$(pwd)/build/"
rm -rf $LARCHBUILDPATH
mkdir $LARCHBUILDPATH
export LD_LIBRARY_PATH="$LARCHBUILDPATH/tbb_cmake_build/tbb_cmake_build_subdir_debug"

cd $LARCHBUILDPATH

cmake ..
# /matsen/cmake-3.27.3-linux-x86_64/bin/cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j16 larch-test
make -j16 larch-usher
make -j16 merge
