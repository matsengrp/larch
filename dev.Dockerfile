FROM ubuntu:22.04
# FROM continuumio/miniconda3:latest

# OPTIONS
ARG NUM_THREADS="4"
RUN export NUM_THREADS=${NUM_THREADS}

RUN apt -y update \
  && apt -y upgrade

# Install required/recommended programs
RUN apt -y install --no-install-recommends \
  wget \
  ssh \
  git \
  vim \
  nano \
  perl \
  black \
  clang-format \
  clang-tidy \
  less
RUN apt -y install --no-install-recommends \
  cmake \
  protobuf-compiler \
  automake \
  autoconf \
  libtool \
  nasm \
  yasm
RUN apt -y install --no-install-recommends \
  wget \
  git \
  ca-certificates \
  make \
  g++ \
  mpi-default-dev \
  libboost-dev \
  libboost-program-options-dev \
  libboost-filesystem-dev \
  libboost-date-time-dev \
  libboost-iostreams-dev

# Copy repo
WORKDIR /app
COPY . /app

# Install larch
WORKDIR /app
RUN rm -rf /app/build
WORKDIR /app/build
RUN cmake -DCMAKE_BUILD_TYPE=Debug ..
RUN make -j${NUM_THREADS}
RUN cp /app/build/larch-usher /usr/local/bin
RUN cp /app/build/larch-test /usr/local/bin
RUN cp

# Working directory
WORKDIR /data

# Start a bash shell when the container launches
CMD ["/bin/bash"]
