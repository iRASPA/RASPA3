FROM rockylinux/rockylinux:10.0

RUN \
  dnf update -y && \
  dnf install -y 'dnf-command(config-manager)' && \
  dnf config-manager --set-enabled crb && \
  dnf install -y --setopt=install_weak_deps=False  git wget binutils binutils-devel && \
  dnf install -y --setopt=install_weak_deps=False  rpm-build xz && \
  dnf install -y --setopt=install_weak_deps=False  openssl-devel && \
  dnf install -y --setopt=install_weak_deps=False  epel-release && \
  dnf install -y --setopt=install_weak_deps=False  libomp zlib-static && \
  dnf install -y perl && \
  dnf install -y --setopt=install_weak_deps=False  lapack-devel lapack64 blas64 && \
  dnf --enablerepo=devel install -y --setopt=install_weak_deps=False python3 python3-devel python3-pybind11 && \
  dnf --enablerepo=devel install -y --setopt=install_weak_deps=False pybind11-devel && \
  dnf --enablerepo=devel install -y --setopt=install_weak_deps=False fftw-devel && \
  dnf --enablerepo=devel install -y --setopt=install_weak_deps=False ocl-icd-devel && \
  wget https://github.com/llvm/llvm-project/releases/download/llvmorg-22.1.1/LLVM-22.1.1-Linux-ARM64.tar.xz && \
  xz -d -v LLVM-22.1.1-Linux-ARM64.tar.xz && \
  tar -xvf LLVM-22.1.1-Linux-ARM64.tar -C /usr --strip-components=1 && \
  rm -f LLVM-22.1.1-Linux-ARM64.tar && \
  ldconfig && \
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
