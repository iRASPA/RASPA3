FROM ubuntu:24.04

RUN \
  apt-get update && \
  DEBIAN_FRONTEND=noninteractive TZ=Europe/Amsterdam apt-get -y install tzdata && \
  apt-get install -y --no-install-recommends git ca-certificates cmake ninja-build libomp5 binutils-dev && \
  apt-get install -y --no-install-recommends llvm lld clang clang-tools clang-tidy libc++-dev libc++abi-dev libomp-18-dev && \
  apt-get install -y --no-install-recommends python3 pybind11-dev python3-pybind11 && \
  apt-get install -y --no-install-recommends liblapack64-dev libblas64-dev && \
  apt-get install -y --no-install-recommends libhdf5-cpp-103 libhdf5-dev python3-h5py python3-tables && \
  apt-get clean

