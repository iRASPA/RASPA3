FROM opensuse/tumbleweed

RUN \
  zypper --non-interactive dup && \
  zypper --non-interactive update  && \
  zypper install -y --no-recommends wget git ca-certificates cmake ninja libstdc++-devel gcc gcc-c++ binutils-devel && \
  zypper install -y --no-recommends python3 python-pybind11-common-devel python3-pybind11 python3-devel && \
  zypper install -y --no-recommends lapack-devel blas-devel && \
  zypper install -y --no-recommends ocl-icd-devel opencl-headers libOpenCL1 fftw3-devel && \
  zypper remove -y busybox && \
  zypper install -y rpm-build git hdf5-devel && \
  zypper install -y llvm clang clang-tools lld && \
  wget https://github.com/llvm/llvm-project/archive/refs/tags/llvmorg-20.1.2.tar.gz && \
  tar -zxvf llvmorg-20.1.2.tar.gz && \
  cd llvm-project-llvmorg-20.1.2 && \
  cmake -G Ninja -S llvm -B build -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_BUILD_TYPE=Release  -DLLVM_ENABLE_PROJECTS="clang;lld;openmp" -DLLVM_ENABLE_RUNTIMES="compiler-rt;libcxx;libcxxabi;libunwind"  -DLLVM_TARGETS_TO_BUILD=X86 -DLLVM_RUNTIME_TARGETS="x86_64-unknown-linux-gnu" -DLLVM_BINUTILS_INCDIR=/usr/include -DLIBOMP_ENABLE_SHARED=OFF && \
  ninja -C build && \
  ninja -C build install && \
  cd .. && \
  rm -rf llvm* && \
  ldconfig && \
  zypper clean
