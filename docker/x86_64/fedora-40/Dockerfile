FROM fedora:40

RUN \
  dnf update -y && \
  dnf update -y && \
  dnf install -y 'dnf-command(config-manager)' && \
  dnf install -y --setopt=install_weak_deps=False  wget git rpm-build && \
  dnf install -y --setopt=install_weak_deps=False  llvm cmake clang clang-tools-extra ninja-build && \
  dnf install -y --setopt=install_weak_deps=False  libomp-devel lld libomp-devel libcxx libcxxabi libcxx-devel libcxxabi-devel libcxx-static libcxxabi-static && \
  dnf install -y --setopt=install_weak_deps=False  lapack-devel lapack64 blas64 && \
  dnf install -y --setopt=install_weak_deps=False python3 python3-devel python3-pybind11 && \
  dnf install -y --setopt=install_weak_deps=False pybind11-devel && \
  dnf clean all
