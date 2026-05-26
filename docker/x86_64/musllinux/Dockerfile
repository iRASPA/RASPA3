FROM quay.io/pypa/musllinux_1_1_x86_64

RUN \
  apk update && apk upgrade && \
  apk add --no-cache musl-dev build-base && \
  apk add --no-cache opencl-headers opencl-icd-loader-dev libgomp && \
  apk add --no-cache gcc g++ make python3 python3-dev openssl-dev && \
  apk add --no-cache lapack-dev blas-dev openssl-dev wget git binutils-dev && \
  apk add --no-cache opencl-headers opencl-icd-loader fftw-dev && \
  wget https://github.com/openssl/openssl/releases/download/openssl-3.0.8/openssl-3.0.8.tar.gz && \
  tar xvf openssl-3.0.8.tar.gz && \
  cd openssl-3.0*/ && \
  ./config --prefix=/usr --libdir=/usr/lib && \
  make -j 16 && \
  make install && \
  cd .. && \
  wget https://github.com/Kitware/CMake/releases/download/v4.3.2/cmake-4.3.2.tar.gz && \
  tar -zxf cmake-4.3.2.tar.gz && \
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
  wget https://github.com/llvm/llvm-project/archive/refs/tags/llvmorg-22.1.6.tar.gz && \
  tar -zxf llvmorg-22.1.6.tar.gz && \
  cd llvm-project-llvmorg-22.1.6 && \ 
  cmake -G Ninja -S llvm -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=/usr \
  -DLLVM_ENABLE_LLD=OFF \
  -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
  -DLLVM_INCLUDE_TESTS=OFF \
  -DLLVM_ENABLE_PROJECTS="clang;lld" \
  -DLLVM_ENABLE_RUNTIMES="libcxxabi;libcxx;libunwind;compiler-rt" \
  -DLLVM_TARGETS_TO_BUILD="X86" \
  -DLLVM_HOST_TRIPLE=x86_64-unknown-linux-musl \
  -DBUILD_SHARED_LIBS=OFF \
  -DLIBCXX_ENABLE_PIC=ON \
  -DLIBCXXABI_ENABLE_PIC=ON \
  -DLIBUNWIND_ENABLE_PIC=ON \
  -DCMAKE_CXX_FLAGS="-fPIC" \
  -DCMAKE_C_FLAGS="-fPIC" \
  -DLIBCXX_ENABLE_SHARED=OFF -DLIBCXX_ENABLE_STATIC=ON \
  -DLIBCXXABI_ENABLE_SHARED=OFF -DLIBCXXABI_ENABLE_STATIC=ON \
  -DLIBUNWIND_ENABLE_SHARED=OFF -DLIBUNWIND_ENABLE_STATIC=ON \
  -DLIBCXX_HAS_MUSL_LIBC=ON \
  -DLIBCXX_CXX_ABI=libcxxabi \
  -DLIBCXX_USE_COMPILER_RT=ON \
  -DLIBCXX_ENABLE_STATIC_ABI_LIBRARY=ON \
  -DLIBCXXABI_USE_COMPILER_RT=ON \
  -DLIBCXXABI_USE_LLVM_UNWINDER=ON \
  -DLIBCXXABI_ENABLE_STATIC_UNWINDER=ON \
  -DCOMPILER_RT_BUILD_BUILTINS=ON \
  -DCOMPILER_RT_BUILD_SANITIZERS=OFF \
  -DCOMPILER_RT_BUILD_LIBFUZZER=OFF \
  -DCOMPILER_RT_BUILD_XRAY=OFF \
  -DCOMPILER_RT_BUILD_MEMPROF=OFF \
  -DCOMPILER_RT_BUILD_PROFILE=OFF \
  -DCOMPILER_RT_BUILD_ORC=OFF \
  -DCOMPILER_RT_BUILD_CTX_PROFILE=OFF \
  -DCOMPILER_RT_BUILD_GWP_ASAN=OFF \
  -DCOMPILER_RT_DEFAULT_TARGET_ONLY=ON \
  -DLLVM_STATIC_LINK_CXX_STDLIB=ON \
  -DCLANG_DEFAULT_RTLIB=compiler-rt \
  -DCLANG_DEFAULT_CXX_STDLIB=libc++ \
  -DCLANG_DEFAULT_UNWINDLIB=libunwind \
  -DLLVM_BINUTILS_INCDIR=/usr/include && \
  ninja -C build -j 8 && \
  ninja -C build install && \
  cd .. && \
  ln -sf /usr/lib/LLVMgold.so /usr/lib/bfd-plugins/LLVMgold.so 2>/dev/null || \
  cp build/lib/LLVMgold.so /usr/lib/bfd-plugins/ && \
  rm -rf llvm* 
