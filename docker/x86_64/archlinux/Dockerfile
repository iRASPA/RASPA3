FROM archlinux:base

RUN \
  pacman -Sy --noconfirm && \
  pacman -Syu --noconfirm archlinux-keyring && \
  pacman -S --noconfirm vi vim git wget debugedit unzip && \
  pacman -S --noconfirm python3 pybind11 && \
  pacman -S --noconfirm gcc make openmp && \
  pacman -S --noconfirm blas64 lapack64 && \
  pacman -S --noconfirm fftw ocl-icd opencl-headers && \
  pacman -S --noconfirm fakeroot shadow && \
  pacman -S --needed --noconfirm sudo && \
  useradd builduser -m && \
  passwd -d builduser && \
  printf 'builduser ALL=(ALL) ALL\n' | tee -a /etc/sudoers 
