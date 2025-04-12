FROM opensuse/tumbleweed

RUN \
  zypper --non-interactive dup && \
  zypper --non-interactive update  && \
  zypper install -y --no-recommends wget git ca-certificates cmake ninja libstdc++-devel gcc gcc-c++ binutils-devel && \
  zypper install -y --no-recommends python3 python-pybind11-common-devel python3-pybind11 python3-devel && \
  zypper install -y --no-recommends lapack-devel blas-devel && \
  zypper install -y --no-recommends ocl-icd-devel fftw3-devel && \
  zypper remove -y busybox && \
  zypper install -y rpm-build git hdf5-devel && \
  zypper install -y llvm clang clang-tools libc++-devel libc++abi-devel libomp-devel lld && \
  ldconfig && \
  zypper clean
