FROM archlinux:latest

RUN \
  pacman -Sy && \
  pacman -S --noconfirm vi vim git wget && \
  pacman -S --noconfirm cmake ninja python3 pybind11 && \
  pacman -S --noconfirm gcc openmp && \
  pacman -S --noconfirm blas64 lapack64 && \
  pacman -S --noconfirm fakeroot shadow && \
  pacman -S --needed --noconfirm sudo && \
  useradd builduser -m && \
  passwd -d builduser && \
  printf 'builduser ALL=(ALL) ALL\n' | tee -a /etc/sudoers && \
  wget https://github.com/llvm/llvm-project/archive/refs/tags/llvmorg-18.1.8.tar.gz && \
  tar -zxvf llvmorg-18.1.8.tar.gz && \
  cd llvm-project-llvmorg-18.1.8 && \
  cmake -G Ninja -S llvm -B build -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_BUILD_TYPE=Release  -DLLVM_ENABLE_PROJECTS="clang;lld" -DLLVM_ENABLE_RUNTIMES="libcxx;libcxxabi;libunwind"  -DLLVM_TARGETS_TO_BUILD=X86 -DLLVM_RUNTIME_TARGETS="x86_64-unknown-linux-gnu" -DLLVM_BINUTILS_INCDIR=/usr/include && \
  ninja -C build && \
  ninja -C build install && \
  rm -rf build && \
  export CC=clang && \
  export CXX=clang++ && \
  cmake -G Ninja -S llvm -B build -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_BUILD_TYPE=Release -DCMAKE_CC_COMPILER=clang  -DCMAKE_CXX_COMPILER=clang++  -DLLVM_ENABLE_PROJECTS="clang;lld;openmp" -DLLVM_ENABLE_RUNTIMES="libcxx;libcxxabi;libunwind"  -DLLVM_TARGETS_TO_BUILD=X86 -DLLVM_RUNTIME_TARGETS="x86_64-unknown-linux-gnu" -DLLVM_BINUTILS_INCDIR=/usr/include && \
  ninja -C build && \
  ninja -C build install && \
  cd .. && \
  rm -rf llvm* && \
  ldconfig && \
  wget https://github.com/HDFGroup/hdf5/releases/download/hdf5_1.14.4.3/hdf5.tar.gz && \
  tar -zxvf hdf5.tar.gz && \
  cd hdf5-1.14.4-3 && \
  cmake -B build -G Ninja -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_C_COMPILER=clang -DCMAKE_INSTALL_PREFIX=/usr -DHDF5_BUILD_CPP_LIB=ON -DBUILD_SHARED_LIBS=OFF -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -std=c++23 -stdlib=libc++"  -DHDF5_ENABLE_NONSTANDARD_FEATURE_FLOAT16=OFF && \
  ninja -C build && \
  ninja -C build install && \
  cd - && \
  rm -rf hdf5* && \
  ldconfig && \
  rm /var/cache/pacman/pkg/* 
