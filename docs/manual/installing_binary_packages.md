# Installing `RASPA` pre-Compiled Packages
\page installing_binary_packages Installing pre-Compiled Packages

RASPA 3 is a modern Monte‑Carlo simulation package that uses C++23.
Pre‑built binaries for mac and windows are **statically linked** and therefore do not require any external runtime libraries.
If you only want to *use* RASPA 3, download a pre‑built package.
If you want to *develop* or *debug* RASPA 3, follow the build‑from‑source instructions.

All official releases are published on the **GitHub “Releases” tab**.
Choose the file that matches your operating system, CPU architecture and AVX capability.

| Platform        | File suffix | Example                                           |
| --------------- | ----------- | ------------------------------------------------- |
| Conda           | –           | install with `conda install raspa3 raspalib`      |
| Debian / Ubuntu | `.deb`      | `raspa‑{VERSION}‑{OS}‑{ARCH}‑{AVX}.deb`      |
| RPM‑based Linux | `.rpm`      | `raspa‑{VERSION}‑{OS}‑{ARCH}‑{AVX}.rpm` |
| macOS           | `.pkg`      | `raspa‑{VERSION}‑mac‑{ARCH}.pkg`                      |
| Windows         | `.exe`      | `raspa‑{VERSION}‑windows‑{ARCH}‑{AVX}.exe`       |

## Table of Contents
1. [Mac-Apple-Silicon](#Mac-arm64)
2. [Windows-arm64](#Windows-arm64)
3. [Ubuntu-25-arm64](#Ubuntu-25-arm64)
4. [Ubuntu-24-arm64](#Ubuntu-24-arm64)
5. [Debian-13-arm64](#Debian-13-arm64)
6. [Mac-Intel-avx2](#Mac-intel-avx2)
7. [Mac-Intel-avx512](#Mac-intel-avx512)
8. [Windows-Intel-avx2](#Windows-intel-avx2)
9. [Windows-Intel-avx512](#Windows-intel-avx512)
10. [Ubuntu-25-amd64-avx2](#Ubuntu-25-amd64-avx2)
11. [Ubuntu-25-amd64-avx512](#Ubuntu-25-amd64-avx512)
12. [Ubuntu-24-amd64-avx2](#Ubuntu-24-amd64-avx2)
13. [Ubuntu-24-amd64-avx512](#Ubuntu-24-amd64-avx512)
14. [Ubuntu-22-amd64-avx2](#Ubuntu-22-amd64-avx2)
15. [Ubuntu-22-amd64-avx512](#Ubuntu-22-amd64-avx512)
16. [Ubuntu-20-amd64-avx2](#Ubuntu-20-amd64-avx2)
17. [Ubuntu-20-amd64-avx512](#Ubuntu-20-amd64-avx512)
18. [Debian-13-amd64-avx2](#Debian-13-amd64-avx2)
19. [Debian-13-amd64-avx512](#Debian-13-amd64-avx512)
20. [Debian-12-amd64-avx2](#Debian-13-amd64-avx2)
21. [Debian-12-amd64-avx512](#Debian-12-amd64-avx512)
22. [Debian-11-amd64-avx2](#Debian-11-amd64-avx2)
23. [Debian-11-amd64-avx512](#Debian-11-amd64-avx512)
24. [Debian-10-amd64-avx2](#Debian-10-amd64-avx2)
25. [Debian-10-amd64-avx512](#Debian-10-amd64-avx512)
26. [Archlinux-amd64-avx2](#Archlinux-amd64-avx2)
27. [Archlinux-amd64-avx512](#Archlinux-amd64-avx512)
28. [Almalinux-9-amd64-avx2](#Almalinux-9-amd64-avx2)
29. [Almalinux-9-amd64-avx512](#Almalinux-9-amd64-avx512)
30. [Almalinux-8-amd64-avx2](#Almalinux-8-amd64-avx2)
31. [Almalinux-8-amd64-avx512](#Almalinux-8-amd64-avx512)
32. [Redhat-9-amd64-avx2](#Redhat-9-amd64-avx2)
33. [Redhat-9-amd64-avx512](#Redhat-9-amd64-avx512)
34. [Redhat-8-amd64-avx2](#Redhat-8-amd64-avx2)
35. [Redhat-8-amd64-avx512](#Redhat-8-amd64-avx512)
36. [Redhat-7-amd64-avx2](#Redhat-7-amd64-avx2)
37. [Redhat-7-amd64-avx512](#Redhat-7-amd64-avx512)
38. [Redhat-6-amd64-avx2](#Redhat-6-amd64-avx2)
39. [Redhat-6-amd64-avx512](#Redhat-6-amd64-avx512)
40. [Fedora-41-amd64-avx2](#Fedora-41-amd64-avx2)
41. [Fedora-41-amd64-avx512](#Fedora-41-amd64-avx512)
42. [Fedora-40-amd64-avx2](#Fedora-40-amd64-avx2)
43. [Fedora-40-amd64-avx512](#Fedora-40-amd64-avx512)
44. [Fedora-39-amd64-avx2](#Fedora-39-amd64-avx2)
45. [Fedora-39-amd64-avx512](#Fedora-39-amd64-avx512)
46. [Fedora-38-amd64-avx2](#Fedora-38-amd64-avx2)
47. [Fedora-38-amd64-avx512](#Fedora-38-amd64-avx512)
48. [Fedora-37-amd64-avx2](#Fedora-37-amd64-avx2)
49. [Fedora-37-amd64-avx512](#Fedora-37-amd64-avx512)
50. [Fedora-36-amd64-avx2](#Fedora-36-amd64-avx2)
51. [Fedora-36-amd64-avx512](#Fedora-36-amd64-avx512)
52. [Fedora-35-amd64-avx2](#Fedora-35-amd64-avx2)
53. [Fedora-35-amd64-avx512](#Fedora-35-amd64-avx512)
54. [OpenSUSE-Tumbleweed-amd64-avx2](#OpenSUSE-Tumbleweed-amd64-avx2)
55. [OpenSUSE-Tumbleweed-amd64-avx512](#OpenSUSE-Tumbleweed-amd64-avx512)
56. [OpenSUSE-15.6-amd64-avx2](#OpenSUSE-15.6-amd64-avx2)
57. [OpenSUSE-15.6-amd64-avx512](#OpenSUSE-15.6-amd64-avx512)
58. [OpenSUSE-15.5-amd64-avx2](#OpenSUSE-15.5-amd64-avx2)
59. [OpenSUSE-15.5-amd64-avx512](#OpenSUSE-15.5-amd64-avx512)
60. [OpenSUSE-15.4-amd64-avx2](#OpenSUSE-15.4-amd64-avx2)
61. [OpenSUSE-15.4-amd64-avx512](#OpenSUSE-15.4-amd64-avx512)
62. [OpenSUSE-15.3-amd64-avx2](#OpenSUSE-15.3-amd64-avx2)
63. [OpenSUSE-15.3-amd64-avx512](#OpenSUSE-15.3-amd64-avx512)
64. [OpenSUSE-15.2-amd64-avx2](#OpenSUSE-15.2-amd64-avx2)
65. [OpenSUSE-15.2-amd64-avx512](#OpenSUSE-15.2-amd64-avx512)

### ARM64

#### Ubuntu-25-arm64 <a name="Ubuntu-25-arm64"></a>
dependencies: libblas64-3, liblapack64-3 , libgfortran5, libgcc-s1, libomp5, libc6, zlib1g, ocl-icd-libopencl1, libfftw3-double3
tested on: arm64v8/ubuntu:25.04
```
apt update
apt install wget 
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa_3.0.13_arm64_ubuntu-25.deb
apt-get install ./raspa_3.0.13_arm64_ubuntu-25.deb
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Ubuntu-24-arm64 <a name="Ubuntu-25-arm64"></a>
dependencies: libblas64-3, liblapack64-3 , libgfortran5, libgcc-s1, libomp5, libc6, zlib1g, ocl-icd-libopencl1, libfftw3-double3
tested on: arm64v8/ubuntu:24.04
```
apt update
apt install wget 
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa_3.0.13_arm64_ubuntu-24.deb
apt-get install ./raspa_3.0.13_arm64_ubuntu-24.deb
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Debian-13-arm64 <a name="Debian-13-arm64"></a>
dependencies: libblas64-3, liblapack64-3 , libgfortran5, libgcc-s1, libomp5, libc6, zlib1g, ocl-icd-libopencl1, libfftw3-double3
tested on: arm64v8/debian:13
```
apt update
apt install wget 
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa_3.0.13_arm64_debian-13.deb
apt-get install ./raspa_3.0.13_arm64_debian-13.deb
/usr/share/raspa3/tests/unit_tests_raspakit
```

----------------------------------------------------------------------------------

### Intel

#### Ubuntu-25-amd64 with AVX2 support <a name="Ubuntu-25-amd64-avx2"></a>

dependencies:
tested on: ubuntu:25.04
```
apt update
apt install wget 
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa_3.0.13_amd64_ubuntu-25-core-avx2.deb
apt-get install ./raspa_3.0.13_amd64_ubuntu-25-core-avx2.deb
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Ubuntu-25-amd64 with AVX512 support <a name="Ubuntu-25-amd64-avx512"></a>

dependencies: libblas64-3, liblapack64-3 , libquadmath0, libgfortran5, libgcc-s1, libomp5, libc6, zlib1g, ocl-icd-libopencl1, libfftw3-double3
tested on: ubuntu:25.04
```
apt update
apt install wget 
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa_3.0.13_amd64_ubuntu-25_skylake-avx512.deb
apt-get install ./raspa_3.0.13_amd64_ubuntu-25_skylake-avx512.deb
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Ubuntu-24-amd64 with AVX2 support <a name="Ubuntu-24-amd64-avx2"></a>

dependencies: libblas64-3, liblapack64-3 , libquadmath0, libgfortran5, libgcc-s1, libomp5, libc6, zlib1g, ocl-icd-libopencl1, libfftw3-double3
tested on: ubuntu:24.04
```
apt update
apt install wget 
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa_3.0.13_amd64_ubuntu-24-core-avx2.deb
apt-get install ./raspa_3.0.13_amd64_ubuntu-24-core-avx2.deb
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Ubuntu-24-amd64 with AVX512 support <a name="Ubuntu-24-amd64-avx512"></a>

dependencies: libblas64-3, liblapack64-3 , libquadmath0, libgfortran5, libgcc-s1, libomp5, libc6, zlib1g, ocl-icd-libopencl1, libfftw3-double3
tested on: ubuntu:24.04
```
apt update
apt install wget 
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa_3.0.13_amd64_ubuntu-24_skylake-avx512.deb
apt-get install ./raspa_3.0.13_amd64_ubuntu-24_skylake-avx512.deb
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Ubuntu-22-amd64 with AVX2 support <a name="Ubuntu-22-amd64-avx2"></a>

dependencies: libblas64-3, liblapack64-3 , libquadmath0, libgfortran5, libgcc-s1, libomp5, libc6, zlib1g, ocl-icd-libopencl1, libfftw3-double3
tested on: ubuntu:22.04
```
apt update
apt install wget 
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa_3.0.13_amd64_ubuntu-22-core-avx2.deb
apt-get install ./raspa_3.0.13_amd64_ubuntu-22-core-avx2.deb
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Ubuntu-22-amd64 with AVX512 support <a name="Ubuntu-22-amd64-avx512"></a>

dependencies: libblas64-3, liblapack64-3 , libquadmath0, libgfortran5, libgcc-s1, libomp5, libc6, zlib1g, ocl-icd-libopencl1, libfftw3-double3
tested on: ubuntu:22.04
```
apt update
apt install wget 
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa_3.0.13_amd64_ubuntu-22_skylake-avx512.deb
apt-get install ./raspa_3.0.13_amd64_ubuntu-22_skylake-avx512.deb
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Ubuntu-20-amd64 with AVX2 support <a name="Ubuntu-20-amd64-avx2"></a>

dependencies: libblas64-3, liblapack64-3 , libquadmath0, libgfortran5, libgcc-s1, libomp5, libc6, zlib1g, ocl-icd-libopencl1, libfftw3-double3
tested on: ubuntu:20.04
```
apt update
apt install wget 
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa_3.0.13_amd64_ubuntu-20-core-avx2.deb
apt-get install ./raspa_3.0.13_amd64_ubuntu-20-core-avx2.deb
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Ubuntu-20-amd64 with AVX512 support <a name="Ubuntu-20-amd64-avx512"></a>

dependencies: libblas64-3, liblapack64-3 , libquadmath0, libgfortran5, libgcc-s1, libomp5, libc6, zlib1g, ocl-icd-libopencl1, libfftw3-double3
tested on: ubuntu:20.04
```
apt update
apt install wget 
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa_3.0.13_amd64_ubuntu-20_skylake-avx512.deb
apt-get install ./raspa_3.0.13_amd64_ubuntu-20_skylake-avx512.deb
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Debian-13-amd64 with AVX2 support <a name="Debian-13-amd64-avx2"></a>

dependencies: libblas64-3, liblapack64-3 , libquadmath0, libgfortran5, libgcc-s1, libomp5, libc6, zlib1g, ocl-icd-libopencl1, libfftw3-double3
tested on: debian:13
```
apt update
apt install wget 
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa_3.0.13_amd64_debian-13_core-avx2.deb
apt-get install ./raspa_3.0.13_amd64_debian-13_core-avx2.deb
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Debian-13-amd64 with AVX512 support <a name="Debian-13-amd64-avx512"></a>

dependencies: libblas64-3, liblapack64-3 , libquadmath0, libgfortran5, libgcc-s1, libomp5, libc6, zlib1g, ocl-icd-libopencl1, libfftw3-double3
tested on: debian:13
```
apt update
apt install wget 
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa_3.0.13_amd64_debian-13_skylake-avx512.deb
apt-get install ./raspa_3.0.13_amd64_debian-13_skylake-avx512.deb
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Debian-12-amd64 with AVX2 support <a name="Debian-12-amd64-avx2"></a>

dependencies: libblas64-3, liblapack64-3 , libquadmath0, libgfortran5, libgcc-s1, libomp5, libc6, zlib1g, ocl-icd-libopencl1, libfftw3-double3
tested on: debian:12
```
apt update
apt install wget 
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa_3.0.13_amd64_debian-12_core-avx2.deb
apt-get install ./raspa_3.0.13_amd64_debian-12_core-avx2.deb
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Debian-12-amd64 with AVX512 support <a name="Debian-12-amd64-avx512"></a>

dependencies: libblas64-3, liblapack64-3 , libquadmath0, libgfortran5, libgcc-s1, libomp5, libc6, zlib1g, ocl-icd-libopencl1, libfftw3-double3
tested on: debian:12
```
apt update
apt install wget 
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa_3.0.13_amd64_debian-12_skylake-avx512.deb
apt-get install ./raspa_3.0.13_amd64_debian-12_skylake-avx512.deb
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Debian-11-amd64 with AVX2 support <a name="Debian-11-amd64-avx2"></a>

dependencies: libblas64-3, liblapack64-3 , libquadmath0, libgfortran5, libgcc-s1, libomp5, libc6, zlib1g, ocl-icd-libopencl1, libfftw3-double3
tested on: debian:11
```
apt update
apt install wget 
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa_3.0.13_amd64_debian-11_core-avx2.deb
apt-get install ./raspa_3.0.13_amd64_debian-11_core-avx2.deb
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Debian-11-amd64 with AVX512 support <a name="Debian-11-amd64-avx512"></a>

dependencies: libblas64-3, liblapack64-3 , libquadmath0, libgfortran5, libgcc-s1, libomp5, libc6, zlib1g, ocl-icd-libopencl1, libfftw3-double3
tested on: debian:11
```
apt update
apt install wget 
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa_3.0.13_amd64_debian-11_skylake-avx512.deb
apt-get install ./raspa_3.0.13_amd64_debian-11_skylake-avx512.deb
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Debian-10-amd64 with AVX2 support <a name="Debian-10-amd64-avx2"></a>

dependencies: libblas3, liblapack3 , libquadmath0, libgfortran3, libgcc1, libomp5, libc6, zlib1g, ocl-icd-libopencl1, libfftw3-double3
tested on: debian:10
workaround: echo "deb http://archive.debian.org/debian stretch main contrib non-free" > /etc/apt/sources.list
```
apt update
apt install wget 
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa_3.0.13_amd64_debian-10_core-avx2.deb
apt-get install ./raspa_3.0.13_amd64_debian-10_core-avx2.deb
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Debian-10-amd64 with AVX512 support <a name="Debian-10-amd64-avx512"></a>

dependencies: libblas3, liblapack3 , libquadmath0, libgfortran3, libgcc1, libomp5, libc6, zlib1g, ocl-icd-libopencl1, libfftw3-double3
tested on: debian:10
workaround: echo "deb http://archive.debian.org/debian stretch main contrib non-free" > /etc/apt/sources.list
```
apt update
apt install wget 
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa_3.0.13_amd64_debian-10_skylake-avx512.deb
apt-get install ./raspa_3.0.13_amd64_debian-10_skylake-avx512.deb
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Archlinux-amd64 with AVX2 support <a name="Archlinux-amd64-avx2"></a>

dependencies: automatically handled when creating package
tested on: archlinux:latest
```
pacman -Sy
pacman -S wget
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1-x86_64-core-avx2.pkg.tar.zst
pacman -U ./raspa-3.0.13-1-x86_64-core-avx2.pkg.tar.zst
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Archlinux-amd64 with AVX512 support <a name="Archlinux-amd64-avx512"></a>

dependencies: automatically handled when creating package
tested on: archlinux:latest
```
pacman -Sy
pacman -S wget
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1-x86_64-skylake-avx512.pkg.tar.zst
pacman -U ./raspa-3.0.13-1-x86_64-skylake-avx512.pkg.tar.zst
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Almalinux-9-amd64 with AVX2 support <a name="Almalinux-9-amd64-avx2"></a>

dependencies: blas64, lapack64, libgfortran, libquadmath,  libomp, glibc, libgcc, zlib, ocl-icd, fftw3
tested on: almalinux/9-base
```
yum install https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.almalinux.el9.x86_64.core-avx2.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Almalinux-9-amd64 with AVX512 support <a name="Almalinux-9-amd64-avx512"></a>

dependencies: blas64, lapack64, libgfortran, libquadmath,  libomp, glibc, libgcc, zlib, ocl-icd, fftw3
tested on: almalinux/9-base
```
yum install https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.almalinux.el9.x86_64.skylake-avx512.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Almalinux-8-amd64 with AVX2 support <a name="Almalinux-8-amd64-avx2"></a>

dependencies: blas, lapack, libgfortran, libquadmath, libomp, glibc, libgcc, zlib, ocl-icd, fftw3
tested on: almalinux/8-base
```
yum install https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.almalinux.el8.x86_64.core-avx2.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Almalinux-8-amd64 with AVX512 support <a name="Almalinux-8-amd64-avx512"></a>

dependencies: blas, lapack, libgfortran, libquadmath, libomp, glibc, libgcc, zlib, ocl-icd, fftw3
tested on: almalinux/8-base
```
yum install https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.almalinux.el8.x86_64.skylake-avx512.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Redhat-9-amd64 with AVX2 support <a name="Redhat-9-amd64-avx2"></a>

dependencies: blas64, lapack64, libgfortran, libquadmath,  libomp, glibc, libgcc, zlib, ocl-icd, fftw3
tested on: rockylinux/rockylinux:9
```
yum install https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.el9.x86_64.core-avx2.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Redhat-9-amd64 with AVX512 support <a name="Redhat-9-amd64-avx512"></a>

dependencies: blas64, lapack64, libgfortran, libquadmath,  libomp, glibc, libgcc, zlib, ocl-icd, fftw3
tested on: rockylinux/rockylinux:9
```
yum install https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.el9.x86_64.skylake-avx512.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Redhat-8-amd64 with AVX2 support <a name="Redhat-8-amd64-avx2"></a>

dependencies: blas, lapack, libgfortran, libquadmath, libomp, glibc, libgcc, zlib, ocl-icd, fftw3
tested on: rockylinux/rockylinux:8
```
yum install https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.el8.x86_64.core-avx2.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Redhat-8-amd64 with AVX512 support <a name="Redhat-8-amd64-avx512"></a>

dependencies: blas, lapack, libgfortran, libquadmath, libomp, glibc, libgcc, zlib, ocl-icd, fftw3
tested on: rockylinux/rockylinux:8
```
yum install https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.el8.x86_64.skylake-avx512.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Redhat-7-amd64 with AVX2 support <a name="Redhat-7-amd64-avx2"></a>

dependencies: blas, lapack, libgfortran, libquadmath,  glibc, libgcc, zlib, ocl-icd, fftw3
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
yum install https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.el7.x86_64.core-avx2.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Redhat-7-amd64 with AVX512 support <a name="Redhat-7-amd64-avx512"></a>

dependencies: blas, lapack, libgfortran, libquadmath,  glibc, libgcc, zlib, ocl-icd, fftw3
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
yum install https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.el7.x86_64.skylake-avx512.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Redhat-6-amd64 with AVX2 support <a name="Redhat-6-amd64-avx2"></a>

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
yum install --nogpgcheck  https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.el6.x86_64.core-avx2.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Redhat-6-amd64 with AVX512 support <a name="Redhat-6-amd64-avx512"></a>

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
yum install --nogpgcheck https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.el6.x86_64.skylake-avx512.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Fedora-41-amd64 with AVX2 support <a name="Fedora-41-amd64-avx2"></a>

dependencies: blas64, lapack64, libgfortran, libquadmath,  libomp, glibc, libgcc, zlib-ng-compat, fftw-libs-double
tested on: fedora:41
```
dnf install https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.fc41.x86_64.core-avx2.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Fedora-41-amd64 with AVX512 support <a name="Fedora-41-amd64-avx512"></a>

dependencies: blas64, lapack64, libgfortran, libquadmath,  libomp, glibc, libgcc, zlib-ng-compat, fftw-libs-double
tested on: fedora:41
```
dnf install https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.fc41.x86_64.skylake-avx512.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Fedora-40-amd64 with AVX2 support <a name="Fedora-40-amd64-avx2"></a>

dependencies: blas64, lapack64, libgfortran, libquadmath,  libomp, glibc, libgcc, zlib-ng-compat, fftw-libs-double
tested on: fedora:40
```
dnf install https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.fc40.x86_64.core-avx2.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Fedora-40-amd64 with AVX512 support <a name="Fedora-40-amd64-avx512"></a>

dependencies: blas64, lapack64, libgfortran, libquadmath,  libomp, glibc, libgcc, zlib-ng-compat, fftw-libs-double
tested on: fedora:40
```
dnf install https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.fc40.x86_64.skylake-avx512.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Fedora-39-amd64 with AVX2 support <a name="Fedora-39-amd64-avx2"></a>

dependencies: blas64, lapack64, libgfortran, libquadmath,  libomp, glibc, libgcc, zlib, fftw-libs-double
tested on: fedora:39
```
dnf install https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.fc39.x86_64.core-avx2.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Fedora-39-amd64 with AVX512 support <a name="Fedora-39-amd64-avx512"></a>

dependencies: blas64, lapack64, libgfortran, libquadmath,  libomp, glibc, libgcc, zlib, fftw-libs-double
tested on: fedora:39
```
dnf install https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.fc39.x86_64.skylake-avx512.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Fedora-38-amd64 with AVX2 support <a name="Fedora-38-amd64-avx2"></a>

dependencies: blas64, lapack64, libgfortran, libquadmath,  libomp, glibc, libgcc, zlib, fftw-libs-double
tested on: fedora:38
```
dnf install https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.fc38.x86_64.core-avx2.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Fedora-38-amd64 with AVX512 support <a name="Fedora-38-amd64-avx512"></a>

dependencies: blas64, lapack64, libgfortran, libquadmath,  libomp, glibc, libgcc, zlib, fftw-libs-double
tested on: fedora:38
```
dnf install https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.fc38.x86_64.skylake-avx512.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Fedora-37-amd64 with AVX2 support <a name="Fedora-37-amd64-avx2"></a>

dependencies: blas64, lapack64, libgfortran, libquadmath,  libomp, glibc, libgcc, zlib, fftw-libs-double
tested on: fedora:37
```
dnf install https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.fc37.x86_64.core-avx2.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Fedora-37-amd64 with AVX512 support <a name="Fedora-37-amd64-avx512"></a>

dependencies: blas64, lapack64, libgfortran, libquadmath,  libomp, glibc, libgcc, zlib, fftw-libs-double
tested on: fedora:37
```
dnf install https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.fc37.x86_64.skylake-avx512.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Fedora-36-amd64 with AVX2 support <a name="Fedora-36-amd64-avx2"></a>

dependencies: blas64, lapack64, libgfortran, libquadmath,  libomp, glibc, libgcc, zlib, fftw-libs-double
tested on: fedora:36
```
dnf install https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.fc36.x86_64.core-avx2.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Fedora-36-amd64 with AVX512 support <a name="Fedora-36-amd64-avx512"></a>

dependencies: blas64, lapack64, libgfortran, libquadmath,  libomp, glibc, libgcc, zlib, fftw-libs-double
tested on: fedora:36
```
dnf install https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.fc36.x86_64.skylake-avx512.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Fedora-35-amd64 with AVX2 support <a name="Fedora-35-amd64-avx2"></a>

dependencies: blas64, lapack64, libgfortran, libquadmath,  libomp, glibc, libgcc, zlib, fftw-libs-double
tested on: fedora:35
```
dnf install https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.fc35.x86_64.core-avx2.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### Fedora-35-amd64 with AVX512 support <a name="Fedora-35-amd64-avx512"></a>

dependencies: blas64, lapack64, libgfortran, libquadmath,  libomp, glibc, libgcc, zlib, fftw-libs-double
tested on: fedora:35
```
dnf install https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.fc35.x86_64.skylake-avx512.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### OpenSUSE-Tumbleweed-amd64 with AVX2 support <a name="OpenSUSE-Tumbleweed-amd64-avx2"></a>

dependencies: libblas3, liblapack3, libgfortran5, libquadmath0, libLLVM17, libedit0, libomp17-devel, glibc, libgcc_s1, zlib, ocl-icd, fftw3
tested on: opensuse/tumbleweed
```
zypper install https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.opensuse-tumbleweed.x86_64.core-avx2.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### OpenSUSE-Tumbleweed-amd64 with AVX512 support <a name="OpenSUSE-Tumbleweed-amd64-avx512"></a>

dependencies: libblas3, liblapack3, libgfortran5, libquadmath0, libLLVM17, libedit0, libomp17-devel, glibc, libgcc_s1, zlib, ocl-icd, fftw3
tested on: opensuse/tumbleweed
```
zypper install https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.opensuse-tumbleweed.x86_64.skylake-avx512.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### OpenSUSE-15.6-amd64 with AVX2 support <a name="OpenSUSE-15.6-amd64-avx2"></a>

dependencies: libhdf5_cpp103, libhdf5-103, libblas3, liblapack3, libgfortran4, libquadmath0, libLLVM17, libedit0, libomp17-devel, glibc, libgcc_s1, zlib, ocl-icd, fftw3
tested on: opensuse/leap:15.6
```
zypper install https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.opensuse-leap-15.6.x86_64.core-avx2.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### OpenSUSE-15.6-amd64 with AVX512 support <a name="OpenSUSE-15.6-amd64-avx512"></a>

dependencies: libhdf5_cpp103, libhdf5-103, libblas3, liblapack3, libgfortran4, libquadmath0, libLLVM17, libedit0, libomp17-devel, glibc, libgcc_s1, zlib, ocl-icd, fftw3
tested on: opensuse/leap:15.6
```
zypper install https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.opensuse-leap-15.6.x86_64.skylake-avx512.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### OpenSUSE-15.5-amd64 with AVX2 support <a name="OpenSUSE-15.5-amd64-avx2"></a>

dependencies: libblas3, liblapack3, libgfortran4, libquadmath0, libLLVM17, libedit0, libomp17-devel, glibc, libgcc_s1, zlib, ocl-icd, fftw3
tested on: opensuse/leap:15.5
```
zypper install https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.opensuse-leap-15.5.x86_64.core-avx2.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### OpenSUSE-15.5-amd64 with AVX512 support <a name="OpenSUSE-15.5-amd64-avx512"></a>

dependencies: libblas3, liblapack3, libgfortran4, libquadmath0, libLLVM17, libedit0, libomp17-devel, glibc, libgcc_s1, zlib, ocl-icd, fftw3
tested on: opensuse/leap:15.5
```
zypper install https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.opensuse-leap-15.5.x86_64.skylake-avx512.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### OpenSUSE-15.4-amd64 with AVX2 support <a name="OpenSUSE-15.4-amd64-avx2"></a>

dependencies: libblas3, liblapack3, libgfortran4, libquadmath0, libLLVM15, libedit0, libomp15-devel, glibc, libgcc_s1, zlib, ocl-icd, fftw3
tested on: opensuse/leap:15.4
```
zypper install https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.opensuse-leap-15.4.x86_64.core-avx2.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### OpenSUSE-15.4-amd64 with AVX512 support <a name="OpenSUSE-15.4-amd64-avx512"></a>

dependencies: libblas3, liblapack3, libgfortran4, libquadmath0, libLLVM15, libedit0, libomp15-devel, glibc, libgcc_s1, zlib, ocl-icd, fftw3
tested on: opensuse/leap:15.4
```
zypper install https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.opensuse-leap-15.4.x86_64.skylake-avx512.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### OpenSUSE-15.3-amd64 with AVX2 support <a name="OpenSUSE-15.3-amd64-avx2"></a>

dependencies: libblas3, liblapack3, libgfortran4, libquadmath0, libLLVM12, libomp12-devel, glibc, libgcc_s1, zlib, ocl-icd, fftw3
tested on: opensuse/leap:15.3
```
zypper install https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.opensuse-leap-15.3.x86_64.core-avx2.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### OpenSUSE-15.3-amd64 with AVX512 support <a name="OpenSUSE-15.3-amd64-avx512"></a>

dependencies: libblas3, liblapack3, libgfortran4, libquadmath0, libLLVM12, libomp12-devel, glibc, libgcc_s1, zlib, ocl-icd, fftw3
tested on: opensuse/leap:15.3
```
zypper install https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.opensuse-leap-15.3.x86_64.skylake-avx512.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### OpenSUSE-15.2-amd64 with AVX2 support <a name="OpenSUSE-15.2-amd64-avx2"></a>

dependencies: libblas3, liblapack3, libgfortran4, libquadmath0, libLLVM9, libomp9-devel, glibc, libgcc_s1, zlib, ocl-icd, fftw3
tested on: opensuse/leap:15.2
```
zypper install https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.opensuse-leap-15.2.x86_64.core-avx2.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```

#### OpenSUSE-15.2-amd64 with AVX512 support <a name="OpenSUSE-15.2-amd64-avx512"></a>

dependencies: libblas3, liblapack3, libgfortran4, libquadmath0, libLLVM9, libomp9-devel, glibc, libgcc_s1, zlib, ocl-icd, fftw3
tested on: opensuse/leap:15.2
```
zypper install https://github.com/iRASPA/RASPA3/releases/download/v3.0.13/raspa-3.0.13-1.opensuse-leap-15.2.x86_64.skylake-avx512.rpm
/usr/share/raspa3/tests/unit_tests_raspakit
```
