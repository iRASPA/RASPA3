FROM quay.io/pypa/manylinux_2_28_x86_64

RUN \
  dnf update -y && \
  dnf install -y --setopt=install_weak_deps=False clang libomp-devel libomp make python3-lit && \
  dnf install -y --setopt=install_weak_deps=False python3 python3-devel python3.11-pybind11-devel && \
  dnf install -y --setopt=install_weak_deps=False lapack-devel && \
  dnf install -y --setopt=install_weak_deps=False lapack-devel lapack64 blas64  && \
  dnf install -y --setopt=install_weak_deps=False wget git rpm-build && \
  dnf install -y --setopt=install_weak_deps=False binutils-devel && \
  dnf install -y --setopt=install_weak_deps=False fftw-devel && \
  dnf install -y --setopt=install_weak_deps=False perl-Data-Dumper perl-IPC-Cmd perl-Pod-Html && \
  dnf install -y --setopt=install_weak_deps=False ocl-icd-devel && \
  wget https://github.com/openssl/openssl/releases/download/openssl-3.0.8/openssl-3.0.8.tar.gz && \
  tar xvf openssl-3.0.8.tar.gz && \
  cd openssl-3.0*/ && \
  ./config && \
  make -j 16 && \
  make install && \
  cd .. && \
  export PATH=/usr/local/openssl/bin:${PATH}  && \
  export LD_LIBRARY_PATH=/usr/local/openssl/lib:${LD_LIBRARY_PATH} && \
  wget https://github.com/Kitware/CMake/releases/download/v4.3.2/cmake-4.3.2.tar.gz && \
  tar -zxvf cmake-4.3.2.tar.gz && \
  cd cmake-4.3.2 && \
  cmake -S . -B build && \
  cmake --build build && \
  cmake --install build && \
  cd .. &&\
  rm -rf cmake-4.3.2 && \
  rm cmake-4.3.2.tar.gz && \
  wget https://github.com/ninja-build/ninja/archive/refs/tags/v1.13.2.tar.gz && \
  tar -zxvf v1.13.2.tar.gz && \
  cd ninja-1.13.2 && \
  cmake -Bbuild-cmake && \
  cmake --build build-cmake && \
  cmake --install build-cmake && \
  cd .. && \
  rm -rf ninja-1.13.2 && \
  rm -f v1.13.2.tar.gz && \
  wget https://github.com/llvm/llvm-project/archive/refs/tags/llvmorg-22.1.5.tar.gz && \
  tar -zxvf llvmorg-22.1.5.tar.gz && \
  cd llvm-project-llvmorg-22.1.5 && \
  cmake -G Ninja -S llvm -B build -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_BUILD_TYPE=Release  -DLLVM_ENABLE_PROJECTS="clang;lld;openmp" -DLLVM_ENABLE_RUNTIMES="compiler-rt;libcxx;libcxxabi;libunwind"  -DLLVM_TARGETS_TO_BUILD=X86 -DLLVM_RUNTIME_TARGETS="x86_64-unknown-linux-gnu" -DLLVM_BINUTILS_INCDIR=/usr/include && \
  ninja -C build -j4 && \
  ninja -C build install && \
  cd .. && \
  rm -rf llvm* && \
  ldconfig && \
  dnf clean all
