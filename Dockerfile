FROM ubuntu:22.04 AS stage_1
# FROM continuumio/miniconda3:latest
ENV PATH="/miniconda/condabin:${PATH}"
ARG PATH="/miniconda/condabin:${PATH}"

# OPTIONS
ARG CMAKE_NUM_THREADS="4"
RUN export CMAKE_NUM_THREADS=${CMAKE_NUM_THREADS}

RUN apt-get -y update \
  && apt-get -y upgrade

# Install required/recommended programs
RUN apt-get -y install --no-install-recommends \
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
RUN apt-get -y install --no-install-recommends \
  build-essential \
  cmake \
  zlib1g-dev \
  python3 \
  python3-pip \
  python3-pytest
RUN apt-get clean

# Install conda
WORKDIR /tmp
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
  bash miniconda.sh -b -p /miniconda && \
  rm miniconda.sh
ENV PATH /miniconda/condabin:miniconda/bin:$PATH
RUN conda init
SHELL ["conda", "run", "-n", "base", "/bin/bash", "-c"]

# # Install larch-dev environment
# WORKDIR /tmp
# COPY environment.yml .
# RUN conda env create -f environment.yml
# RUN rm environment.yml
# SHELL ["conda", "run", "-n", "larch-dev", "/bin/bash", "-c"]
# RUN echo "conda activate larch-dev"

# # Install larch environment
RUN conda create -n larch
RUN conda install --channel "conda-forge" --update-deps --override-channels \
  cmake \
  make \
  cxx-compiler \
  openmpi \
  openmpi-mpicc \
  openmpi-mpicxx \
  boost-cpp \
  automake \
  autoconf \
  libtool \
  yasm \
  ucx \
  zlib
SHELL ["conda", "run", "-n", "larch", "/bin/bash", "-c"]
RUN echo "conda activate larch"

# Copy repo
WORKDIR /larch-repo
COPY . /larch-repo
RUN rm -rf /larch-repo/build

# Install larch
WORKDIR /larch
RUN cmake -DCMAKE_BUILD_TYPE=Debug /larch-repo
RUN make -j${CMAKE_NUM_THREADS}
# RUN echo "export PATH=/app/build:$PATH" >> ~/.bashrc
RUN ln -s larch-usher /usr/local/bin
RUN ln -s merge /usr/local/bin
RUN ln -s dag2dot /usr/local/bin
RUN ln -s larch-test /usr/local/bin

# Cleanup
RUN rm -rf /larch-repo
RUN rm

# Working directory
WORKDIR /data
