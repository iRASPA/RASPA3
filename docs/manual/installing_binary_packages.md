# Installing `RASPA` pre-Compiled Packages
\page installing_binary_packages Installing pre-Compiled Packages

RASPA 3 is a modern Monte‑Carlo simulation package that uses C++23.
Pre‑built binaries for mac and windows are **statically linked** and therefore do not require any external runtime libraries.
If you only want to *use* RASPA 3, download a pre‑built package.
If you want to *develop* or *debug* RASPA 3, follow the build‑from‑source instructions.

All official releases are published on the **GitHub “Releases” tab**.
Choose the file that matches your operating system and CPU architecture.
There is one package per operating system and architecture: it is built for the baseline instruction set and runs on any processor of that architecture.

| Platform        | File suffix | Example                                           |
| --------------- | ----------- | ------------------------------------------------- |
| Conda           | –           | install with `conda install raspa3 raspalib`      |
| Debian / Ubuntu | `.deb`      | `raspa‑{VERSION}‑{OS}‑{ARCH}.deb`      |
| RPM‑based Linux | `.rpm`      | `raspa‑{VERSION}‑{OS}‑{ARCH}.rpm` |
| macOS           | `.pkg`      | `raspa‑{VERSION}‑mac‑{ARCH}.pkg`                      |
| Windows         | `.exe`      | `raspa‑{VERSION}‑windows‑{ARCH}.exe`       |

## Table of Contents
1. [Ubuntu-25-arm64](#Ubuntu-25-arm64)
2. [Ubuntu-24-arm64](#Ubuntu-24-arm64)
3. [Debian-13-arm64](#Debian-13-arm64)
4. [Ubuntu-25-amd64](#Ubuntu-25-amd64)
5. [Ubuntu-24-amd64](#Ubuntu-24-amd64)
6. [Ubuntu-22-amd64](#Ubuntu-22-amd64)
7. [Ubuntu-20-amd64](#Ubuntu-20-amd64)
8. [Debian-13-amd64](#Debian-13-amd64)
9. [Debian-12-amd64](#Debian-12-amd64)
10. [Debian-11-amd64](#Debian-11-amd64)
11. [Debian-10-amd64](#Debian-10-amd64)
12. [Archlinux-amd64](#Archlinux-amd64)
13. [Almalinux-9-amd64](#Almalinux-9-amd64)
14. [Almalinux-8-amd64](#Almalinux-8-amd64)
15. [Redhat-9-amd64](#Redhat-9-amd64)
16. [Redhat-8-amd64](#Redhat-8-amd64)
17. [Redhat-7-amd64](#Redhat-7-amd64)
18. [Redhat-6-amd64](#Redhat-6-amd64)
19. [Fedora-41-amd64](#Fedora-41-amd64)
20. [Fedora-40-amd64](#Fedora-40-amd64)
21. [Fedora-39-amd64](#Fedora-39-amd64)
22. [Fedora-38-amd64](#Fedora-38-amd64)
23. [Fedora-37-amd64](#Fedora-37-amd64)
24. [Fedora-36-amd64](#Fedora-36-amd64)
25. [Fedora-35-amd64](#Fedora-35-amd64)
26. [OpenSUSE-Tumbleweed-amd64](#OpenSUSE-Tumbleweed-amd64)
27. [OpenSUSE-15.6-amd64](#OpenSUSE-15.6-amd64)
28. [OpenSUSE-15.5-amd64](#OpenSUSE-15.5-amd64)
29. [OpenSUSE-15.4-amd64](#OpenSUSE-15.4-amd64)
30. [OpenSUSE-15.3-amd64](#OpenSUSE-15.3-amd64)
31. [OpenSUSE-15.2-amd64](#OpenSUSE-15.2-amd64)

### ARM64

#### Ubuntu-25-arm64 <a name="Ubuntu-25-arm64"></a>
dependencies: libblas64-3, liblapack64-3, libgfortran5, libgcc-s1, libc6, zlib1g, ocl-icd-libopencl1, libfftw3-double3
tested on: arm64v8/ubuntu:25.04
```
apt update
apt install wget 
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.21/raspa_3.0.21_arm64_ubuntu-25.deb
apt-get install ./raspa_3.0.21_arm64_ubuntu-25.deb
/usr/share/raspa3/tests/unit_tests_structurekit
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Ubuntu-24-arm64 <a name="Ubuntu-24-arm64"></a>
dependencies: libblas64-3, liblapack64-3, libgfortran5, libgcc-s1, libc6, zlib1g, ocl-icd-libopencl1, libfftw3-double3
tested on: arm64v8/ubuntu:24.04
```
apt update
apt install wget 
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.21/raspa_3.0.21_arm64_ubuntu-24.deb
apt-get install ./raspa_3.0.21_arm64_ubuntu-24.deb
/usr/share/raspa3/tests/unit_tests_structurekit
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Debian-13-arm64 <a name="Debian-13-arm64"></a>
dependencies: libblas64-3, liblapack64-3, libgfortran5, libgcc-s1, libc6, zlib1g, ocl-icd-libopencl1, libfftw3-double3
tested on: arm64v8/debian:13
```
apt update
apt install wget 
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.21/raspa_3.0.21_arm64_debian-13.deb
apt-get install ./raspa_3.0.21_arm64_debian-13.deb
/usr/share/raspa3/tests/unit_tests_structurekit
/usr/share/raspa3/tests/unit_tests_raspakit
```

----------------------------------------------------------------------------------

### x86_64

#### Ubuntu-25-amd64 <a name="Ubuntu-25-amd64"></a>

dependencies:
tested on: ubuntu:25.04
```
apt update
apt install wget 
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.21/raspa_3.0.21_amd64_ubuntu-25.deb
apt-get install ./raspa_3.0.21_amd64_ubuntu-25.deb
/usr/share/raspa3/tests/unit_tests_structurekit
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Ubuntu-24-amd64 <a name="Ubuntu-24-amd64"></a>

dependencies: libblas64-3, liblapack64-3, libquadmath0, libgfortran5, libgcc-s1, libc6, zlib1g, ocl-icd-libopencl1, libfftw3-double3
tested on: ubuntu:24.04
```
apt update
apt install wget 
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.21/raspa_3.0.21_amd64_ubuntu-24.deb
apt-get install ./raspa_3.0.21_amd64_ubuntu-24.deb
/usr/share/raspa3/tests/unit_tests_structurekit
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Ubuntu-22-amd64 <a name="Ubuntu-22-amd64"></a>

dependencies: libblas64-3, liblapack64-3, libquadmath0, libgfortran5, libgcc-s1, libc6, zlib1g, ocl-icd-libopencl1, libfftw3-double3
tested on: ubuntu:22.04
```
apt update
apt install wget 
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.21/raspa_3.0.21_amd64_ubuntu-22.deb
apt-get install ./raspa_3.0.21_amd64_ubuntu-22.deb
/usr/share/raspa3/tests/unit_tests_structurekit
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Ubuntu-20-amd64 <a name="Ubuntu-20-amd64"></a>

dependencies: libblas64-3, liblapack64-3, libquadmath0, libgfortran5, libgcc-s1, libc6, zlib1g, ocl-icd-libopencl1, libfftw3-double3
tested on: ubuntu:20.04
```
apt update
apt install wget 
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.21/raspa_3.0.21_amd64_ubuntu-20.deb
apt-get install ./raspa_3.0.21_amd64_ubuntu-20.deb
/usr/share/raspa3/tests/unit_tests_structurekit
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Debian-13-amd64 <a name="Debian-13-amd64"></a>

dependencies: libblas64-3, liblapack64-3, libquadmath0, libgfortran5, libgcc-s1, libc6, zlib1g, ocl-icd-libopencl1, libfftw3-double3
tested on: debian:13
```
apt update
apt install wget 
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.21/raspa_3.0.21_amd64_debian-13.deb
apt-get install ./raspa_3.0.21_amd64_debian-13.deb
/usr/share/raspa3/tests/unit_tests_structurekit
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Debian-12-amd64 <a name="Debian-12-amd64"></a>

dependencies: libblas64-3, liblapack64-3, libquadmath0, libgfortran5, libgcc-s1, libc6, zlib1g, ocl-icd-libopencl1, libfftw3-double3
tested on: debian:12
```
apt update
apt install wget 
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.21/raspa_3.0.21_amd64_debian-12.deb
apt-get install ./raspa_3.0.21_amd64_debian-12.deb
/usr/share/raspa3/tests/unit_tests_structurekit
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Debian-11-amd64 <a name="Debian-11-amd64"></a>

dependencies: libblas64-3, liblapack64-3, libquadmath0, libgfortran5, libgcc-s1, libc6, zlib1g, ocl-icd-libopencl1, libfftw3-double3
tested on: debian:11
```
apt update
apt install wget 
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.21/raspa_3.0.21_amd64_debian-11.deb
apt-get install ./raspa_3.0.21_amd64_debian-11.deb
/usr/share/raspa3/tests/unit_tests_structurekit
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Debian-10-amd64 <a name="Debian-10-amd64"></a>

dependencies: libblas3, liblapack3, libquadmath0, libgfortran3, libgcc1, libc6, zlib1g, ocl-icd-libopencl1, libfftw3-double3
tested on: debian:10
workaround: echo "deb http://archive.debian.org/debian stretch main contrib non-free" > /etc/apt/sources.list
```
apt update
apt install wget 
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.21/raspa_3.0.21_amd64_debian-10.deb
apt-get install ./raspa_3.0.21_amd64_debian-10.deb
/usr/share/raspa3/tests/unit_tests_structurekit
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Archlinux-amd64 <a name="Archlinux-amd64"></a>

dependencies: automatically handled when creating package
tested on: archlinux:latest
```
pacman -Sy
pacman -S wget
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.21/raspa-3.0.21-1-x86_64.pkg.tar.zst
pacman -U ./raspa-3.0.21-1-x86_64.pkg.tar.zst
/usr/share/raspa3/tests/unit_tests_structurekit
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Almalinux-9-amd64 <a name="Almalinux-9-amd64"></a>

dependencies: blas64, lapack64, libgfortran, libquadmath, glibc, libgcc, zlib, ocl-icd, fftw3
tested on: almalinux/9-base
```
yum install https://github.com/iRASPA/RASPA3/releases/download/v3.0.21/raspa-3.0.21-1.almalinux.el9.x86_64.rpm
/usr/share/raspa3/tests/unit_tests_structurekit
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Almalinux-8-amd64 <a name="Almalinux-8-amd64"></a>

dependencies: blas, lapack, libgfortran, libquadmath, glibc, libgcc, zlib, ocl-icd, fftw3
tested on: almalinux/8-base
```
yum install https://github.com/iRASPA/RASPA3/releases/download/v3.0.21/raspa-3.0.21-1.almalinux.el8.x86_64.rpm
/usr/share/raspa3/tests/unit_tests_structurekit
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Redhat-9-amd64 <a name="Redhat-9-amd64"></a>

dependencies: blas64, lapack64, libgfortran, libquadmath, glibc, libgcc, zlib, ocl-icd, fftw3
tested on: rockylinux/rockylinux:9
```
yum install https://github.com/iRASPA/RASPA3/releases/download/v3.0.21/raspa-3.0.21-1.el9.x86_64.rpm
/usr/share/raspa3/tests/unit_tests_structurekit
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Redhat-8-amd64 <a name="Redhat-8-amd64"></a>

dependencies: blas, lapack, libgfortran, libquadmath, glibc, libgcc, zlib, ocl-icd, fftw3
tested on: rockylinux/rockylinux:8
```
yum install https://github.com/iRASPA/RASPA3/releases/download/v3.0.21/raspa-3.0.21-1.el8.x86_64.rpm
/usr/share/raspa3/tests/unit_tests_structurekit
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Redhat-7-amd64 <a name="Redhat-7-amd64"></a>

dependencies: blas, lapack, libgfortran, libquadmath, glibc, libgcc, zlib, ocl-icd, fftw3
tested on: centos:centos7.9.2009
```
sed -i -e "s|mirrorlist=|#mirrorlist=|g" /etc/yum.repos.d/*
sed -i -e "s|#baseurl=|baseurl=|g" /etc/yum.repos.d/*
sed -i -e "s|http://mirror.centos.org|https://vault.centos.org|g" /etc/yum.repos.d/*
yum update -y
yum install epel-release
```

```
yum install blas
yum install https://github.com/iRASPA/RASPA3/releases/download/v3.0.21/raspa-3.0.21-1.el7.x86_64.rpm
/usr/share/raspa3/tests/unit_tests_structurekit
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Redhat-6-amd64 <a name="Redhat-6-amd64"></a>

dependencies: blas, lapack, libgfortran, glibc, zlib, ocl-icd, fftw3
tested on: centos:centos6
```
curl https://www.getpagespeed.com/files/centos6-eol.repo --output /etc/yum.repos.d/CentOS-Base.repo 
curl https://www.getpagespeed.com/files/centos6-epel-eol.repo --output /etc/yum.repos.d/epel.repo
yum -y install centos-release-scl
curl https://www.getpagespeed.com/files/centos6-scl-eol.repo --output /etc/yum.repos.d/CentOS-SCLo-scl.repo
curl https://www.getpagespeed.com/files/centos6-scl-rh-eol.repo --output /etc/yum.repos.d/CentOS-SCLo-scl-rh.repo
yum update -y && \
```
```
yum install --nogpgcheck  https://github.com/iRASPA/RASPA3/releases/download/v3.0.21/raspa-3.0.21-1.el6.x86_64.rpm
/usr/share/raspa3/tests/unit_tests_structurekit
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Fedora-41-amd64 <a name="Fedora-41-amd64"></a>

dependencies: blas64, lapack64, libgfortran, libquadmath, glibc, libgcc, zlib-ng-compat, fftw-libs-double
tested on: fedora:41
```
dnf install https://github.com/iRASPA/RASPA3/releases/download/v3.0.21/raspa-3.0.21-1.fc41.x86_64.rpm
/usr/share/raspa3/tests/unit_tests_structurekit
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Fedora-40-amd64 <a name="Fedora-40-amd64"></a>

dependencies: blas64, lapack64, libgfortran, libquadmath, glibc, libgcc, zlib-ng-compat, fftw-libs-double
tested on: fedora:40
```
dnf install https://github.com/iRASPA/RASPA3/releases/download/v3.0.21/raspa-3.0.21-1.fc40.x86_64.rpm
/usr/share/raspa3/tests/unit_tests_structurekit
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Fedora-39-amd64 <a name="Fedora-39-amd64"></a>

dependencies: blas64, lapack64, libgfortran, libquadmath, glibc, libgcc, zlib, fftw-libs-double
tested on: fedora:39
```
dnf install https://github.com/iRASPA/RASPA3/releases/download/v3.0.21/raspa-3.0.21-1.fc39.x86_64.rpm
/usr/share/raspa3/tests/unit_tests_structurekit
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Fedora-38-amd64 <a name="Fedora-38-amd64"></a>

dependencies: blas64, lapack64, libgfortran, libquadmath, glibc, libgcc, zlib, fftw-libs-double
tested on: fedora:38
```
dnf install https://github.com/iRASPA/RASPA3/releases/download/v3.0.21/raspa-3.0.21-1.fc38.x86_64.rpm
/usr/share/raspa3/tests/unit_tests_structurekit
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Fedora-37-amd64 <a name="Fedora-37-amd64"></a>

dependencies: blas64, lapack64, libgfortran, libquadmath, glibc, libgcc, zlib, fftw-libs-double
tested on: fedora:37
```
dnf install https://github.com/iRASPA/RASPA3/releases/download/v3.0.21/raspa-3.0.21-1.fc37.x86_64.rpm
/usr/share/raspa3/tests/unit_tests_structurekit
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Fedora-36-amd64 <a name="Fedora-36-amd64"></a>

dependencies: blas64, lapack64, libgfortran, libquadmath, glibc, libgcc, zlib, fftw-libs-double
tested on: fedora:36
```
dnf install https://github.com/iRASPA/RASPA3/releases/download/v3.0.21/raspa-3.0.21-1.fc36.x86_64.rpm
/usr/share/raspa3/tests/unit_tests_structurekit
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Fedora-35-amd64 <a name="Fedora-35-amd64"></a>

dependencies: blas64, lapack64, libgfortran, libquadmath, glibc, libgcc, zlib, fftw-libs-double
tested on: fedora:35
```
dnf install https://github.com/iRASPA/RASPA3/releases/download/v3.0.21/raspa-3.0.21-1.fc35.x86_64.rpm
/usr/share/raspa3/tests/unit_tests_structurekit
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### OpenSUSE-Tumbleweed-amd64 <a name="OpenSUSE-Tumbleweed-amd64"></a>

dependencies: libblas3, liblapack3, libgfortran5, libquadmath0, libLLVM17, libedit0, glibc, libgcc_s1, zlib, ocl-icd, fftw3
tested on: opensuse/tumbleweed
```
zypper install https://github.com/iRASPA/RASPA3/releases/download/v3.0.21/raspa-3.0.21-1.opensuse-tumbleweed.x86_64.rpm
/usr/share/raspa3/tests/unit_tests_structurekit
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### OpenSUSE-15.6-amd64 <a name="OpenSUSE-15.6-amd64"></a>

dependencies: libhdf5_cpp103, libhdf5-103, libblas3, liblapack3, libgfortran4, libquadmath0, libLLVM17, libedit0, glibc, libgcc_s1, zlib, ocl-icd, fftw3
tested on: opensuse/leap:15.6
```
zypper install https://github.com/iRASPA/RASPA3/releases/download/v3.0.21/raspa-3.0.21-1.opensuse-leap-15.6.x86_64.rpm
/usr/share/raspa3/tests/unit_tests_structurekit
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### OpenSUSE-15.5-amd64 <a name="OpenSUSE-15.5-amd64"></a>

dependencies: libblas3, liblapack3, libgfortran4, libquadmath0, libLLVM17, libedit0, glibc, libgcc_s1, zlib, ocl-icd, fftw3
tested on: opensuse/leap:15.5
```
zypper install https://github.com/iRASPA/RASPA3/releases/download/v3.0.21/raspa-3.0.21-1.opensuse-leap-15.5.x86_64.rpm
/usr/share/raspa3/tests/unit_tests_structurekit
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### OpenSUSE-15.4-amd64 <a name="OpenSUSE-15.4-amd64"></a>

dependencies: libblas3, liblapack3, libgfortran4, libquadmath0, libLLVM15, libedit0, glibc, libgcc_s1, zlib, ocl-icd, fftw3
tested on: opensuse/leap:15.4
```
zypper install https://github.com/iRASPA/RASPA3/releases/download/v3.0.21/raspa-3.0.21-1.opensuse-leap-15.4.x86_64.rpm
/usr/share/raspa3/tests/unit_tests_structurekit
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### OpenSUSE-15.3-amd64 <a name="OpenSUSE-15.3-amd64"></a>

dependencies: libblas3, liblapack3, libgfortran4, libquadmath0, libLLVM12, glibc, libgcc_s1, zlib, ocl-icd, fftw3
tested on: opensuse/leap:15.3
```
zypper install https://github.com/iRASPA/RASPA3/releases/download/v3.0.21/raspa-3.0.21-1.opensuse-leap-15.3.x86_64.rpm
/usr/share/raspa3/tests/unit_tests_structurekit
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### OpenSUSE-15.2-amd64 <a name="OpenSUSE-15.2-amd64"></a>

dependencies: libblas3, liblapack3, libgfortran4, libquadmath0, libLLVM9, glibc, libgcc_s1, zlib, ocl-icd, fftw3
tested on: opensuse/leap:15.2
```
zypper install https://github.com/iRASPA/RASPA3/releases/download/v3.0.21/raspa-3.0.21-1.opensuse-leap-15.2.x86_64.rpm
/usr/share/raspa3/tests/unit_tests_structurekit
/usr/share/raspa3/tests/unit_tests_raspakit
```
