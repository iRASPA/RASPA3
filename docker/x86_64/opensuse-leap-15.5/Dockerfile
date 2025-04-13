FROM opensuse/leap:15.5

RUN \
  zypper update -y && \
  zypper install -y --no-recommends wget git tar ca-certificates libstdc++-devel gcc gcc-c++ binutils-devel openssl-devel && \
  zypper install -y --no-recommends rpm-build && \
  zypper install -y --no-recommends lapack-devel blas-devel && \
  zypper install -y --no-recommends ocl-icd-devel opencl-headers libOpenCL1 fftw3-devel && \
  wget https://www.python.org/ftp/python/3.8.19/Python-3.8.19.tgz && \
  tar xzf Python-3.8.19.tgz && \
  cd Python-3.8.19 && \
  CFLAGS=-fPIC ./configure && \
  make altinstall && \
  ln -sfn /usr/local/bin/python3.8 /usr/bin/python3.8 && \
  cd .. && \
  rm -rf Python-3.8.19* && \
  wget https://github.com/Kitware/CMake/releases/download/v3.31.6/cmake-3.31.6.tar.gz && \
  tar -zxvf cmake-3.31.6.tar.gz && \
  cd cmake-3.31.6 && \
  ./bootstrap && \
  nice -n20 make -j 32 && \
  make install && \
  cd .. &&\
  rm -rf cmake-3.31.6 && \
  rm cmake-3.31.6.tar.gz && \
  wget https://github.com/ninja-build/ninja/archive/refs/tags/v1.12.1.tar.gz && \
  tar -zxvf v1.12.1.tar.gz && \
  cd ninja-1.12.1 && \
  cmake -Bbuild-cmake && \
  cmake --build build-cmake && \
  cmake --install build-cmake && \
  cd .. && \
  rm -rf ninja-1.12.1 && \
  rm  v1.12.1.tar.gz && \
  wget https://github.com/llvm/llvm-project/archive/refs/tags/llvmorg-20.1.2.tar.gz && \
  tar -zxvf llvmorg-20.1.2.tar.gz && \
  cd llvm-project-llvmorg-20.1.2 && \
  cmake -G Ninja -S llvm -B build -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_BUILD_TYPE=Release  -DLLVM_ENABLE_PROJECTS="clang;lld;openmp" -DLLVM_ENABLE_RUNTIMES="compiler-rt;libcxx;libcxxabi;libunwind"  -DLLVM_TARGETS_TO_BUILD=X86 -DLLVM_RUNTIME_TARGETS="x86_64-unknown-linux-gnu" -DLLVM_BINUTILS_INCDIR=/usr/include && \
  ninja -C build && \
  ninja -C build install && \
  cd .. && \
  rm -rf llvm* && \
  ldconfig && \
  wget https://github.com/HDFGroup/hdf5/releases/download/hdf5_1.14.4.3/hdf5.tar.gz && \
  tar -zxvf hdf5.tar.gz && \
  cd hdf5-1.14.4-3 &&\
  cmake -B build -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_C_COMPILER=clang -DCMAKE_INSTALL_PREFIX=/usr -DHDF5_BUILD_CPP_LIB=ON -DBUILD_SHARED_LIBS=OFF -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -std=c++23 -stdlib=libc++"  -DHDF5_ENABLE_NONSTANDARD_FEATURE_FLOAT16=OFF && \
  cd build && \
  make -j 16 && \
  make install && \
  cd ../.. && \
  rm -rf hdf5-1.14.4-3 && \
  ldconfig && \
  zypper clean
