#!/bin/sh
sed -i '' 's/const target& this_ = \*this;/const class target\& this_ = *this;/' deps/usher/mutation_detailed.pb.cc
