FROM ubuntu:22.04

WORKDIR /home
COPY . .

RUN apt -y update
RUN apt -y upgrade
RUN apt -y install --no-install-recommends wget git ca-certificates make g++ mpi-default-dev libboost-dev libboost-program-options-dev libboost-filesystem-dev libboost-date-time-dev libboost-iostreams-dev libtool yasm

RUN wget https://github.com/Kitware/CMake/releases/download/v3.27.3/cmake-3.27.3-linux-x86_64.tar.gz
RUN tar -xvf /home/cmake-3.27.3-linux-x86_64.tar.gz
RUN git clone https://github.com/matsengrp/larch.git
#RUN git submodule set-url deps/usher https://github.com/matsengrp/usher.git
WORKDIR /home/larch/deps
RUN echo "[submodule \"deps/usher\"]\n\tpath = deps/usher\n\turl = https://github.com/matsengrp/usher.git\n\tbranch = add-usher-callback" > ../.gitmodules
RUN git submodule init
RUN git submodule update
WORKDIR /home/larch/build
RUN mkdir -p /home/larch/build
RUN /home/cmake-3.27.3-linux-x86_64/bin/cmake -DCMAKE_BUILD_TYPE=Debug ..
ENV LD_LIBRARY_PATH /home/larch/build/tbb_cmake_build/tbb_cmake_build_subdir_debug:$LD_LIBRARY_PATH
RUN make larch-usher
RUN make merge
WORKDIR /home
ENV runlarchusher /home/larch/build/larch-usher
RUN mkdir -p data
