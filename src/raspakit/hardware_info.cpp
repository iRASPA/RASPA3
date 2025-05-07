module;

#ifdef USE_LEGACY_HEADERS
#include <chrono>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <print>
#include <random>
#include <sstream>
#include <type_traits>
#endif

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <VersionHelpers.h>
#include <intrin.h>
#include <sysinfoapi.h>
#include <winsock2.h>
#include <ws2tcpip.h>
#pragma comment(lib, "ws2_32.lib")
#endif

#if defined(__linux__) || defined(__linux)
#include <sys/utsname.h>
#endif

#if defined(__APPLE__)
#include <sys/sysctl.h>
#include <sys/types.h>
#endif

module hardware_info;

#ifndef USE_LEGACY_HEADERS
import <iostream>;
import <random>;
import <chrono>;
import <sstream>;
import <type_traits>;
import <print>;
#endif

import stringutils;
import json;

std::string HardwareInfo::writeInfo()
{
  std::ostringstream stream;
  // see what compiler is used
#if defined(__GNUC__)
#define COMPILER_STRING "gcc " __VERSION__
#elif defined(__INTEL_COMPILER)
#define COMPILER_STRING "Intel icc " __INTEL_COMPILER_BUILD_DATE
#elif defined(__PGI)
#define COMPILER_STRING "PGI compiler "
#elif defined(__SUNPRO_C)
#define COMPILER_STRING "Sun compiler "
#elif defined(__sgi)
#define COMPILER_STRING "SGI compiler " _COMPILER_VERSION
#elif defined(__HP_aCC)
#define COMPILER_STRING "HP compiler "
#elif defined(__DECC)
#define COMPILER_STRING "HP-Compaq-Digital compiler " __DECC_VER
#elif defined(__xlC__)
#define COMPILER_STRING "IBM compiler " __IBMC__
#elif defined(_MSC_VER)
#define COMPILER_STRING "Microsoft compiler " + std::to_string(_MSC_VER)
#elif defined(__BORLANDC__)
#define COMPILER_STRING "Borland compiler "
#elif defined(__MWERKS__)
#define COMPILER_STRING "Metrowerks CodeWarrior "
#else
#define COMPILER_STRING "unknown compiler/version"
#endif

#if defined(WIN32) || defined(_WIN32) || defined(_WIN64) || defined(__WIN32__) || defined(__NT__)
#if _WIN64
  std::print(stream, "Compiled as a 64-bits application\n");
  std::print(stream, "Compiler: {}\n", COMPILER_STRING);
  std::print(stream, "{}, {}\n\n", "Compile Date = " __DATE__, "Compile Time = " __TIME__);
#else
  std::print(stream, "Compiled as a 32-bits application\n");
  std::print(stream, "Compiler: {}\n", COMPILER_STRING);
  std::print(stream, "{}, {}\n\n", "Compile Date = " __DATE__, "Compile Time = " __TIME__);
#endif
#else
#if defined(__LP64__) || defined(__64BIT__) || defined(_LP64) || (__WORDSIZE == 64)
  std::print(stream, "Compiled as a 64-bits application\n");
  std::print(stream, "Compiler: {}\n", COMPILER_STRING);
  std::print(stream, "{}, {}\n\n", "Compile Date = " __DATE__, "Compile Time = " __TIME__);
#else
  std::print(stream, "Compiled as a 32-bits application\n");
  std::print(stream, "Compiler: {}\n", COMPILER_STRING);
  std::print(stream, "{}, {}\n\n", "Compile Date = " __DATE__, "Compile Time = " __TIME__);
#endif
#endif

  // Get a local time_point with system_clock::duration precision
  //  FIX
  // auto now = std::chrono::zoned_time{ std::chrono::current_zone(), std::chrono::system_clock::now()
  // }.get_local_time();
  auto now = std::chrono::system_clock::now();

  // Get a local time_point with days precision
  auto dp = std::chrono::floor<std::chrono::days>(now);

  // Convert local days-precision time_point to a local {y, m, d} calendar
  std::chrono::year_month_day ymd{dp};
  std::chrono::weekday weekday{dp};
  std::chrono::hh_mm_ss<std::chrono::minutes> time{std::chrono::duration_cast<std::chrono::minutes>(now - dp)};
  std::print(stream, "{}\n", now);

  std::print(stream, "Simulation started on {}, {} {}\n", weekday, ymd.month(), ymd.day());
  std::print(stream, "The start time was {}\n\n", time);

  // get hostname and cpu-info for linux
#if defined(__linux__) || defined(__linux)
  struct utsname uts;
  uname(&uts);
  std::print(stream, "Hostname:    {}\n", std::string(uts.nodename, strlen(uts.nodename)));
  std::print(stream, "OS type:     {} on {}\n", std::string(uts.sysname, strlen(uts.sysname)),
             std::string(uts.machine, strlen(uts.machine)));
  std::print(stream, "OS release:  {}\n", std::string(uts.release, strlen(uts.release)));
  std::print(stream, "OS version:  {}\n\n", std::string(uts.version, strlen(uts.version)));
#endif

#if defined(__CYGWIN__)
  int mib[2];
  size_t len;
  len = sizeof(buffer);

  struct utsname uts;
  uname(&uts);
  std::print(stream, "Hostname:    {}\n", uts.nodename);
  std::print(stream, "OS type:     {} on {}\n", uts.sysname, uts.machine);
  std::print(stream, "OS release:  {}\n", uts.release);
  std::print(stream, "OS version:  {}\n\n", uts.version);
#endif

  // get hostname and cpu-info for mac osx
#if defined(__APPLE__)
  size_t len;
  char cpudata[128], cpumodel[128], hostname[256];
  char osrelease[128], ostype[128], osversion[128];

  len = sizeof(cpudata);
  sysctlbyname("hw.machine", &cpudata, &len, NULL, 0);
  std::print(stream, "Cpu data:    {}\n", std::string(cpudata, cpudata + len - 1));

  len = sizeof(cpumodel);
  sysctlbyname("hw.model", &cpumodel, &len, NULL, 0);
  std::print(stream, "Cpu Model:   {}\n", std::string(cpumodel, cpumodel + len - 1));

  len = sizeof(hostname);
  sysctlbyname("kern.hostname", &hostname, &len, NULL, 0);
  std::print(stream, "Host name:  {}s\n", std::string(hostname, hostname + len - 1));

  len = sizeof(osrelease);
  sysctlbyname("kern.osrelease", &osrelease, &len, NULL, 0);
  std::print(stream, "OS release:  {}\n", std::string(osrelease, osrelease + len - 1));

  len = sizeof(ostype);
  sysctlbyname("kern.ostype", &ostype, &len, NULL, 0);
  std::print(stream, "OS type:     {}\n", std::string(ostype, ostype + len - 1));

  len = sizeof(osversion);
  sysctlbyname("kern.osversion", &osversion, &len, NULL, 0);
  std::print(stream, "OS version:  {}\n\n\n", std::string(osversion, osversion + len - 1));
#endif

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
  SYSTEM_INFO systemInfo;
  GetSystemInfo(&systemInfo);
  OSVERSIONINFOEXW osvi;
  ZeroMemory(&osvi, sizeof(OSVERSIONINFOEXW));
  osvi.dwOSVersionInfoSize = sizeof(osvi);
  if (GetVersionExW((LPOSVERSIONINFOW)&osvi))
  {
    std::print(stream, "[From GetVersionExW] Operating System Version : {}.{} build {}\n", osvi.dwMajorVersion,
               osvi.dwMinorVersion, osvi.dwBuildNumber);
    std::print(stream, "[From GetVersionExW] Service Pack :  {}.{}\n", osvi.wServicePackMajor, osvi.wServicePackMinor);
  }

  if (osvi.dwMajorVersion >= 6)
  {
    if (!IsWindowsXPOrGreater())
    {
      std::print(stream, "[From Version Helper-Functions] Pre-Windows XP\n");
    }
    else if (!IsWindowsXPSP1OrGreater())
    {
      std::print(stream, "[From Version Helper-Functions] Windows XP\n");
    }
    else if (!IsWindowsXPSP2OrGreater())
    {
      std::print(stream, "[From Version Helper-Functions] Windows XP-SP1\n");
    }
    else if (!IsWindowsXPSP3OrGreater())
    {
      std::print(stream, "[From Version Helper-Functions] Windows XP-SP2\n");
    }
    else if (!IsWindowsVistaOrGreater())
    {
      std::print(stream, "[From Version Helper-Functions] Windows XP-SP3\n");
    }
    else if (!IsWindowsVistaSP1OrGreater())
    {
      std::print(stream, "[From Version Helper-Functions] Windows Vista\n");
    }
    else if (!IsWindowsVistaSP2OrGreater())
    {
      std::print(stream, "[From Version Helper-Functions] Windows Vista-SP1\n");
    }
    else if (!IsWindows7OrGreater())
    {
      std::print(stream, "[From Version Helper-Functions] Windows Vista-SP2\n");
    }
    else if (!IsWindows7SP1OrGreater())
    {
      std::print(stream, "[From Version Helper-Functions] Windows 7\n");
    }
    else if (!IsWindows8OrGreater())
    {
      std::print(stream, "[From Version Helper-Functions] Windows 7-SP1\n");
    }
    else if (!IsWindows8Point1OrGreater())
    {
      std::print(stream, "[From Version Helper-Functions] Windows 8\n");
    }
    else if (!IsWindows10OrGreater())
    {
      std::print(stream, "[From Version Helper-Functions] Windows 8-Point1\n");
    }
    else
    {
      std::print(stream, "[From Version Helper-Functions] Windows 10 Or Greater\n");
    }

    if (IsWindowsServer())
    {
      std::print(stream, "[From Version Helper-Functions] Window Server\n");
    }

    if (osvi.dwMajorVersion == 10)
    {
      if (osvi.dwMinorVersion == 0)
      {
        if (osvi.wProductType == VER_NT_WORKSTATION)
        {
          if (osvi.dwBuildNumber >= 22000)
            std::print(stream, "Microsoft Windows 11 ");
          else
            std::print(stream, "Microsoft Windows 10 ");
        }
        else
        {
          if (osvi.dwBuildNumber >= 20348)
            std::print(stream, "Microsoft Windows Server 2022 ");
          else if (osvi.dwBuildNumber >= 17763)
            std::print(stream, "Microsoft Windows Server 2019 ");
          else
            std::print(stream, "Microsoft Windows Server 2016 ");
        }
      }
    }
    else if (osvi.dwMajorVersion == 6)
    {
      if (osvi.dwMinorVersion == 3)
      {
        if (osvi.wProductType == VER_NT_WORKSTATION)
          std::print(stream, "Microsoft Windows 8.1 ");
        else
          std::print(stream, "Microsoft Windows Server 2012 R2 ");
      }
      else if (osvi.dwMinorVersion == 2)
      {
        if (osvi.wProductType == VER_NT_WORKSTATION)
          std::print(stream, "Microsoft Windows 8 ");
        else
          std::print(stream, "Microsoft Windows Server 2012 ");
      }
      else if (osvi.dwMinorVersion == 1)
      {
        if (osvi.wProductType == VER_NT_WORKSTATION)
          std::print(stream, "Microsoft Windows 7 ");
        else
          std::print(stream, "Microsoft Windows Server 2008 R2 ");
      }
      else if (osvi.dwMinorVersion == 0)
      {
        if (osvi.wProductType == VER_NT_WORKSTATION)
          std::print(stream, "Microsoft Windows Vista ");
        else
          std::print(stream, "Microsoft Windows Server 2008 ");
      }
    }
  }

  DWORD dwType;
  if ((GetProductInfo != nullptr) && GetProductInfo(osvi.dwMajorVersion, osvi.dwMinorVersion, 0, 0, &dwType))
  {
    switch (dwType)
    {
      case PRODUCT_ULTIMATE:
        std::print(stream, "Ultimate Edition\n");
        break;
      case PRODUCT_PROFESSIONAL:
        std::print(stream, "Professional\n");
        break;
      case PRODUCT_HOME_PREMIUM:
        std::print(stream, "Home Premium Edition\n");
        break;
      case PRODUCT_HOME_BASIC:
        std::print(stream, "Home Basic Edition\n");
        break;
      case PRODUCT_ENTERPRISE:
        std::print(stream, "Enterprise Edition\n");
        break;
      case PRODUCT_BUSINESS:
        std::print(stream, "Business Edition\n");
        break;
      case PRODUCT_STARTER:
        std::print(stream, "Starter Edition\n");
        break;
      case PRODUCT_CLUSTER_SERVER:
        std::print(stream, "Cluster Server Edition\n");
        break;
      case PRODUCT_DATACENTER_SERVER:
        std::print(stream, "Datacenter Edition\n");
        break;
      case PRODUCT_DATACENTER_SERVER_CORE:
        std::print(stream, "Datacenter Edition (core installation)\n");
        break;
      case PRODUCT_ENTERPRISE_SERVER:
        std::print(stream, "Enterprise Edition\n");
        break;
      case PRODUCT_ENTERPRISE_SERVER_CORE:
        std::print(stream, "Enterprise Edition (core installation)\n");
        break;
      case PRODUCT_ENTERPRISE_SERVER_IA64:
        std::print(stream, "Enterprise Edition for Itanium-based Systems\n");
        break;
      case PRODUCT_SMALLBUSINESS_SERVER:
        std::print(stream, "Small Business Server\n");
        break;
      case PRODUCT_SMALLBUSINESS_SERVER_PREMIUM:
        std::print(stream, "Small Business Server Premium Edition\n");
        break;
      case PRODUCT_STANDARD_SERVER:
        std::print(stream, "Standard Edition\n");
        break;
      case PRODUCT_STANDARD_SERVER_CORE:
        std::print(stream, "Standard Edition (core installation)\n");
        break;
      case PRODUCT_WEB_SERVER:
        std::print(stream, "Web Server Edition\n");
        break;
    }
  }

  WSADATA wsaData;
  if (WSAStartup(MAKEWORD(2, 2), &wsaData) == 0)
  {
    char hostname[128] = "";
    gethostname(hostname, sizeof(hostname));
    struct hostent* ent = gethostbyname(hostname);
    if (ent)
    {
      struct in_addr ip_addr = *(struct in_addr*)(ent->h_addr);
      std::print(stream, "Hostname: {}, ip address: {}\n", hostname, inet_ntoa(ip_addr));
    }
  }

  if (systemInfo.wProcessorArchitecture == PROCESSOR_ARCHITECTURE_INTEL)
    std::print(stream, "Processor Architecture 32-bit\n");
  else if (systemInfo.wProcessorArchitecture == PROCESSOR_ARCHITECTURE_AMD64)
    std::print(stream, "Processor Architecture 64-bit\n");
  else if (systemInfo.wProcessorArchitecture == PROCESSOR_ARCHITECTURE_IA64)
    std::print(stream, "Processor Architecture Intel Itanium\n");
  else if (systemInfo.wProcessorArchitecture == PROCESSOR_ARCHITECTURE_ARM)
    std::print(stream, "Processor Architecture ARM\n");
#if !defined(__MINGW32__) && !defined(__MINGW64__)
  else if (systemInfo.wProcessorArchitecture == PROCESSOR_ARCHITECTURE_ARM64)
    std::print(stream, "Processor Architecture ARM64\n");
#endif

#if defined(__x86_64__) || (_M_64__)
  int CPUInfo[4] = {-1};
  unsigned nExIds, i = 0;
  char CPUBrandString[0x40];
  // Get the information associated with each extended ID.
  __cpuid(CPUInfo, 0x80000000);
  nExIds = CPUInfo[0];
  for (i = 0x80000000; i <= nExIds; ++i)
  {
    __cpuid(CPUInfo, i);
    // Interpret CPU brand string
    if (i == 0x80000002)
      memcpy(CPUBrandString, CPUInfo, sizeof(CPUInfo));
    else if (i == 0x80000003)
      memcpy(CPUBrandString + 16, CPUInfo, sizeof(CPUInfo));
    else if (i == 0x80000004)
      memcpy(CPUBrandString + 32, CPUInfo, sizeof(CPUInfo));
  }
  // string includes manufacturer, model and clockspeed
  std::print(stream, "CPU Type:{}\n", CPUBrandString);
#endif

  MEMORYSTATUSEX statex;
  statex.dwLength = sizeof(statex);
  GlobalMemoryStatusEx(&statex);
  std::print(stream, "Total System Memory : {} MB\n", (statex.ullTotalPhys / 1024) / 1024);

  std::print(stream, "Number of Cores:  {}\n\n", systemInfo.dwNumberOfProcessors);
#endif

  return stream.str();
}

nlohmann::json HardwareInfo::jsonInfo()
{
  nlohmann::json info;
  // see what compiler is used
#if defined(__GNUC__)
#define COMPILER_STRING "gcc " __VERSION__
#elif defined(__INTEL_COMPILER)
#define COMPILER_STRING "Intel icc " __INTEL_COMPILER_BUILD_DATE
#elif defined(__PGI)
#define COMPILER_STRING "PGI compiler "
#elif defined(__SUNPRO_C)
#define COMPILER_STRING "Sun compiler "
#elif defined(__sgi)
#define COMPILER_STRING "SGI compiler " _COMPILER_VERSION
#elif defined(__HP_aCC)
#define COMPILER_STRING "HP compiler "
#elif defined(__DECC)
#define COMPILER_STRING "HP-Compaq-Digital compiler " __DECC_VER
#elif defined(__xlC__)
#define COMPILER_STRING "IBM compiler " __IBMC__
#elif defined(_MSC_VER)
#define COMPILER_STRING "Microsoft compiler " + std::to_string(_MSC_VER)
#elif defined(__BORLANDC__)
#define COMPILER_STRING "Borland compiler "
#elif defined(__MWERKS__)
#define COMPILER_STRING "Metrowerks CodeWarrior "
#else
#define COMPILER_STRING "unknown compiler/version"
#endif

  info["compiler"] = std::format("{}", COMPILER_STRING);
  info["compile_date"] = std::format("{}", __DATE__);
  info["compile_time"] = std::format("{}", __TIME__);

  // Get a local time_point with system_clock::duration precision
  //  FIX
  // auto now = std::chrono::zoned_time{ std::chrono::current_zone(), std::chrono::system_clock::now()
  // }.get_local_time();
  auto now = std::chrono::system_clock::now();

  // Get a local time_point with days precision
  auto dp = std::chrono::floor<std::chrono::days>(now);

  // Convert local days-precision time_point to a local {y, m, d} calendar
  std::chrono::year_month_day ymd{dp};
  std::chrono::weekday weekday{dp};
  std::chrono::hh_mm_ss<std::chrono::minutes> time{std::chrono::duration_cast<std::chrono::minutes>(now - dp)};

  info["start_date"] = std::format("{}, {} {}", weekday, ymd.month(), ymd.day());
  info["start_time"] = std::format("{}", time);

  // get hostname and cpu-info for linux
#if defined(__linux__) || defined(__linux)
  struct utsname uts;
  uname(&uts);

  info["hostname"] = std::string(uts.nodename, strlen(uts.nodename));
  info["os_type"] = std::string(uts.sysname, strlen(uts.sysname));
  info["os_arch"] = std::string(uts.machine, strlen(uts.machine));
  info["os_release"] = std::string(uts.release, strlen(uts.release));
  info["os_version"] = std::string(uts.version, strlen(uts.version));
#endif

#if defined(__CYGWIN__)
  int mib[2];
  size_t len;
  len = sizeof(buffer);

  struct utsname uts;
  uname(&uts);
  info["hostname"] = uts.nodename;
  info["os_type"] = uts.sysname;
  info["os_arch"] = uts.machine;
  info["os_release"] = uts.release;
  info["os_version"] = uts.version;
#endif

  // get hostname and cpu-info for mac osx
#if defined(__APPLE__)
  size_t len;
  char cpudata[128], cpumodel[128], hostname[256];
  char osrelease[128], ostype[128], osversion[128];

  len = sizeof(cpudata);
  sysctlbyname("hw.machine", &cpudata, &len, NULL, 0);
  info["os_arch"] = std::string(cpudata, cpudata + len - 1);

  len = sizeof(cpumodel);
  sysctlbyname("hw.model", &cpumodel, &len, NULL, 0);
  info["cpu_model"] = std::string(cpumodel, cpumodel + len - 1);

  len = sizeof(hostname);
  sysctlbyname("kern.hostname", &hostname, &len, NULL, 0);
  info["hostname"] = std::string(hostname, hostname + len - 1);

  len = sizeof(osrelease);
  sysctlbyname("kern.osrelease", &osrelease, &len, NULL, 0);
  info["os_release"] = std::string(osrelease, osrelease + len - 1);

  len = sizeof(ostype);
  sysctlbyname("kern.ostype", &ostype, &len, NULL, 0);
  info["os_type"] = std::string(ostype, ostype + len - 1);

  len = sizeof(osversion);
  sysctlbyname("kern.osversion", &osversion, &len, NULL, 0);
  info["os_version"] = std::string(osversion, osversion + len - 1);
#endif

  /*
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
  SYSTEM_INFO systemInfo;
  GetSystemInfo(&systemInfo);
  OSVERSIONINFOEXW osvi;
  ZeroMemory(&osvi, sizeof(OSVERSIONINFOEXW));
  osvi.dwOSVersionInfoSize = sizeof(osvi);
  if (GetVersionExW((LPOSVERSIONINFOW)&osvi))
  {
    std::print(stream, "[From GetVersionExW] Operating System Version : {}.{} build {}\n", osvi.dwMajorVersion,
               osvi.dwMinorVersion, osvi.dwBuildNumber);
    std::print(stream, "[From GetVersionExW] Service Pack :  {}.{}\n", osvi.wServicePackMajor,
osvi.wServicePackMinor);
  }

  if (osvi.dwMajorVersion >= 6)
  {
    if (!IsWindowsXPOrGreater())
    {
      std::print(stream, "[From Version Helper-Functions] Pre-Windows XP\n");
    }
    else if (!IsWindowsXPSP1OrGreater())
    {
      std::print(stream, "[From Version Helper-Functions] Windows XP\n");
    }
    else if (!IsWindowsXPSP2OrGreater())
    {
      std::print(stream, "[From Version Helper-Functions] Windows XP-SP1\n");
    }
    else if (!IsWindowsXPSP3OrGreater())
    {
      std::print(stream, "[From Version Helper-Functions] Windows XP-SP2\n");
    }
    else if (!IsWindowsVistaOrGreater())
    {
      std::print(stream, "[From Version Helper-Functions] Windows XP-SP3\n");
    }
    else if (!IsWindowsVistaSP1OrGreater())
    {
      std::print(stream, "[From Version Helper-Functions] Windows Vista\n");
    }
    else if (!IsWindowsVistaSP2OrGreater())
    {
      std::print(stream, "[From Version Helper-Functions] Windows Vista-SP1\n");
    }
    else if (!IsWindows7OrGreater())
    {
      std::print(stream, "[From Version Helper-Functions] Windows Vista-SP2\n");
    }
    else if (!IsWindows7SP1OrGreater())
    {
      std::print(stream, "[From Version Helper-Functions] Windows 7\n");
    }
    else if (!IsWindows8OrGreater())
    {
      std::print(stream, "[From Version Helper-Functions] Windows 7-SP1\n");
    }
    else if (!IsWindows8Point1OrGreater())
    {
      std::print(stream, "[From Version Helper-Functions] Windows 8\n");
    }
    else if (!IsWindows10OrGreater())
    {
      std::print(stream, "[From Version Helper-Functions] Windows 8-Point1\n");
    }
    else
    {
      std::print(stream, "[From Version Helper-Functions] Windows 10 Or Greater\n");
    }

    if (IsWindowsServer())
    {
      std::print(stream, "[From Version Helper-Functions] Window Server\n");
    }

    if (osvi.dwMajorVersion == 10)
    {
      if (osvi.dwMinorVersion == 0)
      {
        if (osvi.wProductType == VER_NT_WORKSTATION)
        {
          if (osvi.dwBuildNumber >= 22000)
            std::print(stream, "Microsoft Windows 11 ");
          else
            std::print(stream, "Microsoft Windows 10 ");
        }
        else
        {
          if (osvi.dwBuildNumber >= 20348)
            std::print(stream, "Microsoft Windows Server 2022 ");
          else if (osvi.dwBuildNumber >= 17763)
            std::print(stream, "Microsoft Windows Server 2019 ");
          else
            std::print(stream, "Microsoft Windows Server 2016 ");
        }
      }
    }
    else if (osvi.dwMajorVersion == 6)
    {
      if (osvi.dwMinorVersion == 3)
      {
        if (osvi.wProductType == VER_NT_WORKSTATION)
          std::print(stream, "Microsoft Windows 8.1 ");
        else
          std::print(stream, "Microsoft Windows Server 2012 R2 ");
      }
      else if (osvi.dwMinorVersion == 2)
      {
        if (osvi.wProductType == VER_NT_WORKSTATION)
          std::print(stream, "Microsoft Windows 8 ");
        else
          std::print(stream, "Microsoft Windows Server 2012 ");
      }
      else if (osvi.dwMinorVersion == 1)
      {
        if (osvi.wProductType == VER_NT_WORKSTATION)
          std::print(stream, "Microsoft Windows 7 ");
        else
          std::print(stream, "Microsoft Windows Server 2008 R2 ");
      }
      else if (osvi.dwMinorVersion == 0)
      {
        if (osvi.wProductType == VER_NT_WORKSTATION)
          std::print(stream, "Microsoft Windows Vista ");
        else
          std::print(stream, "Microsoft Windows Server 2008 ");
      }
    }
  }

  DWORD dwType;
  if ((GetProductInfo != nullptr) && GetProductInfo(osvi.dwMajorVersion, osvi.dwMinorVersion, 0, 0, &dwType))
  {
    switch (dwType)
    {
      case PRODUCT_ULTIMATE:
        std::print(stream, "Ultimate Edition\n");
        break;
      case PRODUCT_PROFESSIONAL:
        std::print(stream, "Professional\n");
        break;
      case PRODUCT_HOME_PREMIUM:
        std::print(stream, "Home Premium Edition\n");
        break;
      case PRODUCT_HOME_BASIC:
        std::print(stream, "Home Basic Edition\n");
        break;
      case PRODUCT_ENTERPRISE:
        std::print(stream, "Enterprise Edition\n");
        break;
      case PRODUCT_BUSINESS:
        std::print(stream, "Business Edition\n");
        break;
      case PRODUCT_STARTER:
        std::print(stream, "Starter Edition\n");
        break;
      case PRODUCT_CLUSTER_SERVER:
        std::print(stream, "Cluster Server Edition\n");
        break;
      case PRODUCT_DATACENTER_SERVER:
        std::print(stream, "Datacenter Edition\n");
        break;
      case PRODUCT_DATACENTER_SERVER_CORE:
        std::print(stream, "Datacenter Edition (core installation)\n");
        break;
      case PRODUCT_ENTERPRISE_SERVER:
        std::print(stream, "Enterprise Edition\n");
        break;
      case PRODUCT_ENTERPRISE_SERVER_CORE:
        std::print(stream, "Enterprise Edition (core installation)\n");
        break;
      case PRODUCT_ENTERPRISE_SERVER_IA64:
        std::print(stream, "Enterprise Edition for Itanium-based Systems\n");
        break;
      case PRODUCT_SMALLBUSINESS_SERVER:
        std::print(stream, "Small Business Server\n");
        break;
      case PRODUCT_SMALLBUSINESS_SERVER_PREMIUM:
        std::print(stream, "Small Business Server Premium Edition\n");
        break;
      case PRODUCT_STANDARD_SERVER:
        std::print(stream, "Standard Edition\n");
        break;
      case PRODUCT_STANDARD_SERVER_CORE:
        std::print(stream, "Standard Edition (core installation)\n");
        break;
      case PRODUCT_WEB_SERVER:
        std::print(stream, "Web Server Edition\n");
        break;
    }
  }

  WSADATA wsaData;
  if (WSAStartup(MAKEWORD(2, 2), &wsaData) == 0)
  {
    char hostname[128] = "";
    gethostname(hostname, sizeof(hostname));
    struct hostent* ent = gethostbyname(hostname);
    if (ent)
    {
      struct in_addr ip_addr = *(struct in_addr*)(ent->h_addr);
      std::print(stream, "Hostname: {}, ip address: {}\n", hostname, inet_ntoa(ip_addr));
    }
  }

  if (systemInfo.wProcessorArchitecture == PROCESSOR_ARCHITECTURE_INTEL)
    std::print(stream, "Processor Architecture 32-bit\n");
  else if (systemInfo.wProcessorArchitecture == PROCESSOR_ARCHITECTURE_AMD64)
    std::print(stream, "Processor Architecture 64-bit\n");
  else if (systemInfo.wProcessorArchitecture == PROCESSOR_ARCHITECTURE_IA64)
    std::print(stream, "Processor Architecture Intel Itanium\n");
  else if (systemInfo.wProcessorArchitecture == PROCESSOR_ARCHITECTURE_ARM)
    std::print(stream, "Processor Architecture ARM\n");
#if !defined(__MINGW32__) && !defined(__MINGW64__)
  else if (systemInfo.wProcessorArchitecture == PROCESSOR_ARCHITECTURE_ARM64)
    std::print(stream, "Processor Architecture ARM64\n");
#endif

#if defined(__x86_64__) || (_M_64__)
  int CPUInfo[4] = {-1};
  unsigned nExIds, i = 0;
  char CPUBrandString[0x40];
  // Get the information associated with each extended ID.
  __cpuid(CPUInfo, 0x80000000);
  nExIds = CPUInfo[0];
  for (i = 0x80000000; i <= nExIds; ++i)
  {
    __cpuid(CPUInfo, i);
    // Interpret CPU brand string
    if (i == 0x80000002)
      memcpy(CPUBrandString, CPUInfo, sizeof(CPUInfo));
    else if (i == 0x80000003)
      memcpy(CPUBrandString + 16, CPUInfo, sizeof(CPUInfo));
    else if (i == 0x80000004)
      memcpy(CPUBrandString + 32, CPUInfo, sizeof(CPUInfo));
  }
  // string includes manufacturer, model and clockspeed
  std::print(stream, "CPU Type:{}\n", CPUBrandString);
#endif

  MEMORYSTATUSEX statex;
  statex.dwLength = sizeof(statex);
  GlobalMemoryStatusEx(&statex);
  std::print(stream, "Total System Memory : {} MB\n", (statex.ullTotalPhys / 1024) / 1024);

  std::print(stream, "Number of Cores:  {}\n\n", systemInfo.dwNumberOfProcessors);
#endif

  return stream.str();

  */
  return info;
}
