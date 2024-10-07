FROM fedora:38

RUN \
  dnf update -y && \
  dnf update -y && \
  dnf install -y 'dnf-command(config-manager)' && \
  dnf install -y --setopt=install_weak_deps=False  wget git rpm-build openssl-devel binutils-devel perl && \
  dnf install -y --setopt=install_weak_deps=False  llvm clang clang-tools-extra ninja-build && \
  dnf install -y --setopt=install_weak_deps=False  lapack-devel lapack64 blas64 && \
  dnf install -y --setopt=install_weak_deps=False python3 python3-devel python3-pybind11 && \
  dnf install -y --setopt=install_weak_deps=False pybind11-devel && \
  wget https://github.com/Kitware/CMake/releases/download/v3.29.3/cmake-3.29.3.tar.gz && \
  tar -zxvf cmake-3.29.3.tar.gz && \
  cd cmake-3.29.3 && \
  ./bootstrap && \
  nice -n20 make -j 32 && \
  make install && \
  cd .. &&\
  rm -rf cmake-3.29.3 && \
  rm cmake-3.29.3.tar.gz && \
  wget https://github.com/llvm/llvm-project/archive/refs/tags/llvmorg-18.1.8.tar.gz && \
  tar -zxvf llvmorg-18.1.8.tar.gz && \
  cd llvm-project-llvmorg-18.1.8 && \
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
  dnf clean all
