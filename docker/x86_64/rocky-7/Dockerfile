FROM centos:centos7.9.2009

RUN \
  sed -i -e 's/mirrorlist=/#mirrorlist=/g' /etc/yum.repos.d/CentOS-* && \
  sed -i -e 's/mirrorlist=/#mirrorlist=/g' /etc/yum.conf && \
  sed -E -i -e 's/#baseurl=http:\/\/mirror.centos.org\/centos\/\$releasever\/([[:alnum:]_-]*)\/\$basearch\//baseurl=https:\/\/vault.centos.org\/7.9.2009\/\1\/\$basearch\//g' /etc/yum.repos.d/CentOS-* && \
  sed -E -i -e 's/#baseurl=http:\/\/mirror.centos.org\/centos\/\$releasever\/([[:alnum:]_-]*)\/\$basearch\//baseurl=https:\/\/vault.centos.org\/7.9.2009\/\1\/\$basearch\//g' /etc/yum.conf && \
  yum update -y && \
  yum install -y https://cmcc.mirrors.ustc.edu.cn/centos-vault/centos/7.9.2009/sclo/x86_64/rh/Packages/d/devtoolset-8-gcc-8.3.1-3.2.el7.x86_64.rpm https://cmcc.mirrors.ustc.edu.cn/centos-vault/centos/7.9.2009/sclo/x86_64/rh/Packages/d/devtoolset-8-gcc-c%2B%2B-8.3.1-3.2.el7.x86_64.rpm https://cmcc.mirrors.ustc.edu.cn/centos-vault/centos/7.9.2009/sclo/x86_64/rh/Packages/d/devtoolset-8-binutils-2.30-55.el7.2.x86_64.rpm https://cmcc.mirrors.ustc.edu.cn/centos-vault/centos/7.9.2009/sclo/x86_64/rh/Packages/d/devtoolset-8-libstdc%2B%2B-devel-8.3.1-3.2.el7.x86_64.rpm https://cmcc.mirrors.ustc.edu.cn/centos-vault/centos/7.9.2009/sclo/x86_64/rh/Packages/d/devtoolset-8-runtime-8.1-1.el7.x86_64.rpm && \
  yum install -y --setopt=install_weak_deps=False lapack-devel lapack64 blas64  && \
  yum install -y --setopt=install_weak_deps=False openssl-devel && \
  yum install -y --setopt=install_weak_deps=False wget git && \
  yum install -y --setopt=install_weak_deps=False rpm-build && \
  yum install -y --setopt=install_weak_deps=False make binutils-devel gcc-c++ libstdc++-devel && \
  yum install -y --setopt=install_weak_deps=False epel-release && \
  yum install -y --setopt=install_weak_deps=False python3 python3-devel && \
  yum -y install perl-Data-Dumper && \
  source /opt/rh/devtoolset-8/enable && \
  wget https://github.com/Kitware/CMake/releases/download/v3.29.3/cmake-3.29.3.tar.gz && \
  tar -zxvf cmake-3.29.3.tar.gz && \
  cd cmake-3.29.3 && \
  ./bootstrap && \
  nice -n20 make -j 32 && \
  make install && \
  cd .. &&\
  rm -rf cmake-3.29.3 && \
  rm cmake-3.29.3.tar.gz && \
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
  wget https://github.com/llvm/llvm-project/archive/refs/tags/llvmorg-18.1.7.tar.gz && \
  tar -zxvf llvmorg-18.1.7.tar.gz && \
  cd llvm-project-llvmorg-18.1.7 && \
  cmake -G Ninja -S llvm -B build -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_BUILD_TYPE=Release -DLLVM_ENABLE_PROJECTS="clang;lld;openmp" -DLLVM_ENABLE_RUNTIMES="libcxx;libcxxabi;libunwind"  -DLLVM_TARGETS_TO_BUILD=X86 -DLLVM_RUNTIME_TARGETS="x86_64-unknown-linux-gnu" -DLLVM_BINUTILS_INCDIR=/usr/include  -DLIBOMP_ENABLE_SHARED=OFF && \
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
  yum clean all
