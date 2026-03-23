FROM ghcr.io/iraspa/debian-10-llvm-20:latest

RUN \
  export CC=clang && \
  export CXX=clang++ && \
  wget https://github.com/openssl/openssl/releases/download/openssl-3.0.8/openssl-3.0.8.tar.gz && \
  tar xvf openssl-3.0.8.tar.gz && \
  cd openssl-3.0*/ && \
  ./config && \
  make -j 16 && \
  make install && \
  cd .. && \
  export PATH=/usr/local/openssl/bin:${PATH}  && \
  export LD_LIBRARY_PATH=/usr/local/openssl/lib:${LD_LIBRARY_PATH} && \
  wget https://github.com/Kitware/CMake/releases/download/v4.3.0/cmake-4.3.0.tar.gz && \
  tar -zxvf cmake-4.3.0.tar.gz && \
  cd cmake-4.3.0 && \
  ./bootstrap && \
  nice -n20 make -j 32 && \
  make install && \
  cd .. &&\
  rm -rf cmake-4.3.0 && \
  rm cmake-4.3.0.tar.gz && \
  wget https://github.com/ninja-build/ninja/archive/refs/tags/v1.13.2.tar.gz && \
  tar -zxvf v1.13.2.tar.gz && \
  cd ninja-1.13.2 && \
  cmake -Bbuild-cmake && \
  cmake --build build-cmake && \
  cmake --install build-cmake && \
  cd .. && \
  rm -rf ninja-1.13.2 && \
  rm  v1.13.2.tar.gz && \
  wget https://github.com/llvm/llvm-project/archive/refs/tags/llvmorg-22.1.1.tar.gz && \
  tar -zxvf llvmorg-22.1.1.tar.gz && \
  cd llvm-project-llvmorg-22.1.1 && \
  cmake -G Ninja -S llvm -B build -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_BUILD_TYPE=Release  -DLLVM_ENABLE_PROJECTS="clang;lld;openmp" -DLLVM_ENABLE_RUNTIMES="libcxx;libcxxabi;libunwind"  -DLLVM_TARGETS_TO_BUILD=X86 -DLLVM_RUNTIME_TARGETS="x86_64-unknown-linux-gnu" -DLLVM_BINUTILS_INCDIR=/usr/include && \
  ninja -C build && \
  ninja -C build install && \
  cd .. && \
  rm -rf llvm* && \
  ldconfig && \
  apt-get clean

