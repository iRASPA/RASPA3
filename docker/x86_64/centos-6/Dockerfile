FROM centos:6

RUN \
  curl https://www.getpagespeed.com/files/centos6-eol.repo --output /etc/yum.repos.d/CentOS-Base.repo && \
  curl https://www.getpagespeed.com/files/centos6-epel-eol.repo --output /etc/yum.repos.d/epel.repo && \
  yum -y install centos-release-scl && \
  curl https://www.getpagespeed.com/files/centos6-scl-eol.repo --output /etc/yum.repos.d/CentOS-SCLo-scl.repo && \
  curl https://www.getpagespeed.com/files/centos6-scl-rh-eol.repo --output /etc/yum.repos.d/CentOS-SCLo-scl-rh.repo && \
  yum update -y && \
  yum install -y --nogpgcheck lapack-devel lapack blas  && \
  yum install -y --nogpgcheck openssl-devel rpm-build && \
  yum install -y --nogpgcheck wget git && \
  yum install -y --nogpgcheck binutils-devel && \
  yum install -y --nogpgcheck gcc openssl-devel bzip2-devel && \
  yum install -y --nogpgcheck libstdc++-devel && \
  yum install -y --nogpgcheck ocl-icd-devel && \
  yum install -y --nogpgcheck fftw-devel && \
  yum install -y --nogpgcheck http://buildlogs-seed.centos.org/centos/6/sclo/x86_64/rh/Packages/d/devtoolset-8-gcc-8.3.1-3.el6.x86_64.rpm http://buildlogs-seed.centos.org/centos/6/sclo/x86_64/rh/Packages/d/devtoolset-8-gcc-c%2B%2B-8.3.1-3.2.el6.x86_64.rpm http://buildlogs-seed.centos.org/centos/6/sclo/x86_64/rh/Packages/d/devtoolset-8-binutils-2.30-55.el6.2.x86_64.rpm http://buildlogs-seed.centos.org/centos/6/sclo/x86_64/rh/Packages/d/devtoolset-8-binutils-devel-2.30-55.el6.2.x86_64.rpm http://buildlogs-seed.centos.org/centos/6/sclo/x86_64/rh/Packages/d/devtoolset-8-runtime-8.1-1.el6.x86_64.rpm http://buildlogs-seed.centos.org/centos/6/sclo/x86_64/rh/Packages/d/devtoolset-8-libstdc%2B%2B-devel-8.3.1-3.2.el6.x86_64.rpm http://buildlogs-seed.centos.org/centos/6/sclo/x86_64/rh/Packages/d/devtoolset-8-libatomic-devel-8.3.1-3.2.el6.x86_64.rpm && \
  source /opt/rh/devtoolset-8/enable && \
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
  rm -f cmake-3.31.6.tar.gz && \
  wget https://github.com/ninja-build/ninja/archive/refs/tags/v1.12.1.tar.gz && \
  tar -zxvf v1.12.1.tar.gz && \
  cd ninja-1.12.1 && \
  cmake -Bbuild-cmake && \
  cmake --build build-cmake && \
  cmake --install build-cmake && \
  cd .. && \
  rm -rf ninja-1.12.1 && \
  rm -f v1.12.1.tar.gz && \
  yum clean all && \
  wget https://github.com/llvm/llvm-project/archive/refs/tags/llvmorg-20.1.2.tar.gz && \
  tar -zxvf llvmorg-20.1.2.tar.gz && \
  cd llvm-project-llvmorg-20.1.2 && \
  cmake -G Ninja -S llvm -B build -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_BUILD_TYPE=Release -DCMAKE_CC_COMPILER=gcc  -DCMAKE_CXX_COMPILER=g++  -DLLVM_ENABLE_PROJECTS="clang;lld" -DLLVM_ENABLE_RUNTIMES="libcxx;libcxxabi;libunwind"  -DLLVM_TARGETS_TO_BUILD=X86 -DLLVM_RUNTIME_TARGETS="x86_64-unknown-linux-gnu" -DLLVM_BINUTILS_INCDIR=/opt/rh/devtoolset-8/root/usr/include/ && \
  ninja -C build || true && \
  sed -i '42d;43d' build/include/c++/v1/__memory/aligned_alloc.h && \
  sed -i -e 's/-lpthread/-lpthread -lrt/g' build/runtimes/runtimes-x86_64-unknown-linux-gnu-bins/build.ninja && \
  ninja -C build && \
  ninja -C build install && \
  cp -af /usr/lib/x86_64-unknown-linux-gnu/libc++* /usr/lib64/ && \
  cp /usr/lib/x86_64-unknown-linux-gnu/libunwind.* /usr/lib64/ && \
  rm -rf build && \
  cmake -G Ninja -S llvm -B build -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_C_COMPILER=clang -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_BUILD_TYPE=Release -DCMAKE_CC_COMPILER=clang  -DCMAKE_CXX_COMPILER=clang++  -DLLVM_ENABLE_PROJECTS="clang;lld;openmp" -DLLVM_ENABLE_RUNTIMES="libcxx;libcxxabi;libunwind"  -DLLVM_TARGETS_TO_BUILD=X86 -DLLVM_RUNTIME_TARGETS="x86_64-unknown-linux-gnu" -DLLVM_BINUTILS_INCDIR=/opt/rh/devtoolset-8/root/usr/include/ -DLIBOMP_ENABLE_SHARED=OFF -DCMAKE_CXX_FLAGS="--stdlib=libc++ -fno-diagnostics-color" -DHAVE_GETPAGESIZE=1 -DCMAKE_EXE_LINKER_FLAGS="-fuse-ld=lld -stdlib=libc++" -DCMAKE_SHARED_LINKER_FLAGS="-fuse-ld=lld -stdlib=libc++" -DCMAKE_EXE_LINKER_FLAGS_INIT="-fuse-ld=lld"  -DCMAKE_MODULE_LINKER_FLAGS_INIT="-fuse-ld=lld" -DCMAKE_SHARED_LINKER_FLAGS_INIT="-fuse-ld=lld" && \
  ninja -C build || true && \
  sed -i '42d;43d' build/include/c++/v1/__memory/aligned_alloc.h && \
  sed -i -e 's/-lpthread/-lpthread -lrt/g' build/runtimes/runtimes-x86_64-unknown-linux-gnu-bins/build.ninja && \
  ninja -C build && \
  rm -rf /usr/lib64/libc++* && \
  rm -rf /usr/lib64/libunwind* && \
  ninja -C build install && \
  cd .. && \
  /bin/rm -rf llvm-* && \
  echo "/usr/lib/x86_64-unknown-linux-gnu/" >> /etc/ld.so.conf.d/libomp.conf && \
  ldconfig && \
  wget https://github.com/HDFGroup/hdf5/releases/download/hdf5_1.14.4.3/hdf5.tar.gz && \
  tar -zxvf hdf5.tar.gz && \
  cd hdf5-1.14.4-3 &&\
  cmake -B build -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_C_COMPILER=clang -DCMAKE_INSTALL_PREFIX=/usr -DHDF5_BUILD_CPP_LIB=ON -DBUILD_SHARED_LIBS=OFF -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -std=c++23 -stdlib=libc++"  -DHDF5_ENABLE_NONSTANDARD_FEATURE_FLOAT16=OFF -DCMAKE_EXE_LINKER_FLAGS_INIT="-fuse-ld=lld"  -DCMAKE_MODULE_LINKER_FLAGS_INIT="-fuse-ld=lld" -DCMAKE_SHARED_LINKER_FLAGS_INIT="-fuse-ld=lld" && \
  cd build && \
  make -j 16 && \
  make install && \
  cd ../.. && \
  rm -rf hdf5-1.14.4-3 && \
  ldconfig && \
  yum clean all 
