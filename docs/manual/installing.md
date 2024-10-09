# Installing `RASPA`
\page installing Installing

Download one of the precompiled packages from

> <https://github.com/iRASPA/RASPA3/releases>

In Figure ![1.1](figure_binary_packages.jpg) one can find the list of packages.
These include packages for `macOS` (both `intel` and `apple silicon`),
windows (both `intel` en `arm64`) and many `Linux` distributions.

![List of available binary packages: `Linux` distributions based on
`Red Hat` and `OpenSUSE` (`rpm`-packages), `Debian` (`deb`-packages),
and `Arch` `Linux` (`tar.zst` packages), installers for `intel`- and
`apple-silicon` macs, and installers for `intel`- and `arm64-Windows`
computers.](introduction/raspa3_releases.png){#fig: binary_packages
width="95%"}

## Archlinux
```
pacman -Sy
pacman -S wget
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.0/raspa-3.0.0-1-x86_64.pkg.tar.zst
pacman -U ./raspa-3.0.0-1-x86_64.pkg.tar.zst
```

## Redhat
- redhat-9: `yum install https://github.com/iRASPA/RASPA3/releases/download/v3.0.0/raspa-3.0.0-1.el9.x86_64.rpm`
- redhat-8: `yum install https://github.com/iRASPA/RASPA3/releases/download/v3.0.0/raspa-3.0.0-1.el8.x86_64.rpm`
- redhat-7: `yum install https://github.com/iRASPA/RASPA3/releases/download/v3.0.0/raspa-3.0.0-1.el7.x86_64.rpm`
- redhat-6: `yum install https://github.com/iRASPA/RASPA3/releases/download/v3.0.0/raspa-3.0.0-1.el6.x86_64.rpm`

## Fedora 
- fedora-40: `dnf install https://github.com/iRASPA/RASPA3/releases/download/v3.0.0/raspa-3.0.0-1.fc40.x86_64.rpm`
- fedora-39: `dnf install https://github.com/iRASPA/RASPA3/releases/download/v3.0.0/raspa-3.0.0-1.fc39.x86_64.rpm`
- fedora-38: `dnf install https://github.com/iRASPA/RASPA3/releases/download/v3.0.0/raspa-3.0.0-1.fc38.x86_64.rpm`
- fedora-37: `dnf install https://github.com/iRASPA/RASPA3/releases/download/v3.0.0/raspa-3.0.0-1.fc37.x86_64.rpm`
- fedora-36: `dnf install https://github.com/iRASPA/RASPA3/releases/download/v3.0.0/raspa-3.0.0-1.fc36.x86_64.rpm`
- fedora-35: `dnf install https://github.com/iRASPA/RASPA3/releases/download/v3.0.0/raspa-3.0.0-1.fc35.x86_64.rpm`


## OpenSUSE Tumbleweed (Signature verification failed: choose 'ignore')
```
  zypper update
  zypper install https://github.com/iRASPA/RASPA3/releases/download/v3.0.0/raspa-3.0.0-1.opensuse-tumbleweed.x86_64.rpm
```

## OpenSUSE Leap (Signature verification failed: choose 'ignore')
- 15.5: `zypper install https://github.com/iRASPA/RASPA3/releases/download/v3.0.0/raspa-3.0.0-1.opensuse-leap-15.5.x86_64.rpm`
- 15.4: `zypper install https://github.com/iRASPA/RASPA3/releases/download/v3.0.0/raspa-3.0.0-1.opensuse-leap-15.4.x86_64.rpm`
- 15.3: `zypper install https://github.com/iRASPA/RASPA3/releases/download/v3.0.0/raspa-3.0.0-1.opensuse-leap-15.3.x86_64.rpm`
- 15.2: `zypper install https://github.com/iRASPA/RASPA3/releases/download/v3.0.0/raspa-3.0.0-1.opensuse-leap-15.2.x86_64.rpm`

## Ubuntu
- ubuntu-24
```
apt update
apt install wget
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.0/raspa_3.0.0_amd64-ubuntu-24.deb
apt install ./raspa_3.0.0_amd64-ubuntu-24.deb
```
- ubuntu-22
```
apt update
apt install wget
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.0/raspa_3.0.0_amd64-ubuntu-22.deb
apt install ./raspa_3.0.0_amd64-ubuntu-22.deb
```
- ubuntu-20
```
apt update
apt install wget
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.0/raspa_3.0.0_amd64-ubuntu-20.deb
apt install ./raspa_3.0.0_amd64-ubuntu-20.deb
```

## Debian
- 12
```
apt update
apt install wget
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.0/raspa_3.0.0_amd64-debian-12.deb
apt install ./raspa_3.0.0_amd64-debian-12.deb
```
- 11
```
apt update
apt install wget
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.0/raspa_3.0.0_amd64-debian-11.deb
apt install ./raspa_3.0.0_amd64-debian-11.deb
```
- 10
```
apt update
apt install wget
wget https://github.com/iRASPA/RASPA3/releases/download/v3.0.0/raspa_3.0.0_amd64-debian-10.deb
apt install ./raspa_3.0.0_amd64-debian-10.deb
```

