module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <complex>
#include <cstddef>
#include <functional>
#include <iostream>
#include <istream>
#include <map>
#include <ostream>
#include <ranges>
#include <string>
#include <tuple>
#include <vector>
#endif

export module archive;

#ifndef USE_LEGACY_HEADERS
import <string>;
import <tuple>;
import <vector>;
import <array>;
import <map>;
import <ranges>;
import <istream>;
import <ostream>;
import <iostream>;
import <algorithm>;
import <bit>;
import <chrono>;
import <complex>;
import <functional>;
#endif

// on linux uint64_t is unsigned long        8
//          size_t   is unsigned long        8
//          size_t is an alias for uint64_t
// on mac   uint64_t is unsigned long long   8
//          size_t   is unsigned long        8
//          size_t is an alias for unsigned long

export template <class STREAM>
class Archive
{
 public:
  Archive(STREAM& stream) : stream(stream) {};

  // a function that can serialize any enum
  template <typename Enum, typename = typename std::enable_if<std::is_enum<Enum>::value>::type>
  Archive& operator<<(const Enum& e)
  {
    *this << static_cast<int64_t>(e);
    return *this;
  }

  // a function that can deserialize any enum from an archive
  template <typename Enum, typename = typename std::enable_if<std::is_enum<Enum>::value>::type>
  Archive& operator>>(Enum& e)
  {
    int64_t v;
    *this >> v;
    e = static_cast<Enum>(v);
    return *this;
  }
  Archive& operator>>(bool& v)
  {
    uint8_t w;
    *this >> w;
    v = static_cast<bool>(w);
    return *this;
  }

  Archive& operator<<(const bool& v)
  {
    *this << static_cast<uint8_t>(v);
    return *this;
  }

  Archive& operator>>(int8_t& v)
  {
    stream.read(std::bit_cast<char*>(&v), sizeof(int8_t));
    if (!stream)
    {
      throw std::runtime_error("malformed data");
    }
    return *this;
  }

  Archive& operator<<(const int8_t& v)
  {
    stream.write(std::bit_cast<const char*>(&v), sizeof(int8_t));
    return *this;
  }

  Archive& operator>>(uint8_t& v)
  {
    stream.read(std::bit_cast<char*>(&v), sizeof(uint8_t));
    if (!stream)
    {
      throw std::runtime_error("malformed data");
    }
    return *this;
  }

  Archive& operator<<(const uint8_t& v)
  {
    stream.write(std::bit_cast<const char*>(&v), sizeof(uint8_t));
    return *this;
  }

  Archive& operator>>(int16_t& v)
  {
    stream.read(std::bit_cast<char*>(&v), sizeof(int16_t));
    if (!stream)
    {
      throw std::runtime_error("malformed data");
    }
    if constexpr (std::endian::native == std::endian::little)
    {
      v = std::byteswap(v);
    }
    return *this;
  }

  Archive& operator<<(const int16_t& v)
  {
    int16_t w{v};
    if constexpr (std::endian::native == std::endian::little)
    {
      w = std::byteswap(w);
    }
    stream.write(std::bit_cast<const char*>(&w), sizeof(int16_t));
    return *this;
  }

  Archive& operator>>(uint16_t& v)
  {
    stream.read(std::bit_cast<char*>(&v), sizeof(uint16_t));
    if (!stream)
    {
      throw std::runtime_error("malformed data");
    }
    if constexpr (std::endian::native == std::endian::little)
    {
      v = std::byteswap(v);
    }
    return *this;
  }

  Archive& operator<<(const uint16_t& v)
  {
    uint16_t w{v};
    if constexpr (std::endian::native == std::endian::little)
    {
      w = std::byteswap(w);
    }
    stream.write(std::bit_cast<const char*>(&w), sizeof(uint16_t));
    return *this;
  }

  Archive& operator>>(int32_t& v)
  {
    stream.read(std::bit_cast<char*>(&v), sizeof(int32_t));
    if (!stream)
    {
      throw std::runtime_error("malformed data");
    }
    if constexpr (std::endian::native == std::endian::little)
    {
      v = std::byteswap(v);
    }
    return *this;
  }

  Archive& operator<<(const int32_t& v)
  {
    int32_t w{v};
    if constexpr (std::endian::native == std::endian::little)
    {
      w = std::byteswap(w);
    }
    stream.write(std::bit_cast<const char*>(&w), sizeof(int32_t));
    return *this;
  }

  Archive& operator>>(uint32_t& v)
  {
    stream.read(std::bit_cast<char*>(&v), sizeof(uint32_t));
    if (!stream)
    {
      throw std::runtime_error("malformed data");
    }
    if constexpr (std::endian::native == std::endian::little)
    {
      v = std::byteswap(v);
    }
    return *this;
  }

  Archive& operator<<(const uint32_t& v)
  {
    uint32_t w{v};
    if constexpr (std::endian::native == std::endian::little)
    {
      w = std::byteswap(w);
    }
    stream.write(std::bit_cast<const char*>(&w), sizeof(uint32_t));
    return *this;
  }

  Archive& operator>>(int64_t& v)
  {
    stream.read(std::bit_cast<char*>(&v), sizeof(int64_t));
    if (!stream)
    {
      throw std::runtime_error("malformed data");
    }
    if constexpr (std::endian::native == std::endian::little)
    {
      v = std::byteswap(v);
    }
    return *this;
  }

  Archive& operator<<(const int64_t& v)
  {
    int64_t w{v};
    if constexpr (std::endian::native == std::endian::little)
    {
      w = std::byteswap(w);
    }
    stream.write(std::bit_cast<const char*>(&w), sizeof(int64_t));
    return *this;
  }

#if defined(__APPLE__) && defined(__MACH__)
  Archive& operator>>(uint64_t& v)
  {
    stream.read(std::bit_cast<char*>(&v), sizeof(uint64_t));
    if (!stream)
    {
      throw std::runtime_error("malformed data");
    }
    if constexpr (std::endian::native == std::endian::little)
    {
      v = std::byteswap(v);
    }
    return *this;
  }

  Archive& operator<<(const uint64_t& v)
  {
    uint64_t w{v};
    if constexpr (std::endian::native == std::endian::little)
    {
      w = std::byteswap(w);
    }
    stream.write(std::bit_cast<const char*>(&w), sizeof(uint64_t));
    return *this;
  }
#endif

  Archive& operator>>(size_t& v)
  {
    stream.read(std::bit_cast<char*>(&v), sizeof(size_t));
    if (!stream)
    {
      throw std::runtime_error("malformed data");
    }
    if constexpr (std::endian::native == std::endian::little)
    {
      v = std::byteswap(v);
    }
    return *this;
  }

  Archive& operator<<(const size_t& v)
  {
    size_t w{v};
    if constexpr (std::endian::native == std::endian::little)
    {
      w = std::byteswap(w);
    }
    stream.write(std::bit_cast<const char*>(&w), sizeof(size_t));
    return *this;
  }

#if defined(__APPLE__) && defined(__MACH__)
  Archive& operator>>(std::make_signed_t<std::size_t>& v)
  {
    stream.read(std::bit_cast<char*>(&v), sizeof(std::make_signed_t<std::size_t>));
    if (!stream)
    {
      throw std::runtime_error("malformed data");
    }
    if constexpr (std::endian::native == std::endian::little)
    {
      v = std::byteswap(v);
    }
    return *this;
  }

  Archive& operator<<(const std::make_signed_t<std::size_t>& v)
  {
    std::make_signed_t<std::size_t> w{v};
    if constexpr (std::endian::native == std::endian::little)
    {
      w = std::byteswap(w);
    }
    stream.write(std::bit_cast<const char*>(&w), sizeof(std::make_signed_t<std::size_t>));
    return *this;
  }
#endif

  Archive& operator>>(double& v)
  {
    stream.read(std::bit_cast<char*>(&v), sizeof(double));
    if (!stream)
    {
      throw std::runtime_error("malformed data");
    }
    if constexpr (std::endian::native == std::endian::little)
    {
      auto value_representation = std::bit_cast<std::array<std::byte, sizeof(double)>>(v);
      std::ranges::reverse(value_representation);
      v = std::bit_cast<double>(value_representation);
    }
    return *this;
  }

  Archive<std::ofstream>& operator<<(const double& v)
  {
    double w{v};
    if constexpr (std::endian::native == std::endian::little)
    {
      auto value_representation = std::bit_cast<std::array<std::byte, sizeof(double)>>(w);
      std::ranges::reverse(value_representation);
      w = std::bit_cast<double>(value_representation);
    }
    stream.write(std::bit_cast<const char*>(&w), sizeof(double));
    return *this;
  }

  Archive& operator>>(std::chrono::duration<double>& v)
  {
    double count{0.0};
    stream.read(std::bit_cast<char*>(&count), sizeof(double));
    if (!stream)
    {
      throw std::runtime_error("malformed data");
    }
    if constexpr (std::endian::native == std::endian::little)
    {
      auto value_representation = std::bit_cast<std::array<std::byte, sizeof(double)>>(count);
      std::ranges::reverse(value_representation);
      count = std::bit_cast<double>(value_representation);
    }
    v = std::chrono::duration<double>(count);
    return *this;
  }

  Archive<std::ofstream>& operator<<(const std::chrono::duration<double>& v)
  {
    double w{v.count()};
    if constexpr (std::endian::native == std::endian::little)
    {
      auto value_representation = std::bit_cast<std::array<std::byte, sizeof(double)>>(w);
      std::ranges::reverse(value_representation);
      w = std::bit_cast<double>(value_representation);
    }
    stream.write(std::bit_cast<const char*>(&w), sizeof(double));
    return *this;
  }

  Archive& operator>>(std::complex<double>& v)
  {
    double real, imag;
    *this >> real;
    *this >> imag;
    v = std::complex<double>(real, imag);
    return *this;
  }

  Archive<std::ofstream>& operator<<(const std::complex<double>& v)
  {
    *this << v.real();
    *this << v.imag();
    return *this;
  }

  template <class T>
  Archive& operator>>(std::optional<T>& v)
  {
    v.reset();
    bool has_value;
    *this >> has_value;
    if (has_value)
    {
      T element;
      *this >> element;
      v = std::move(element);
    }
    return *this;
  }

  template <class T>
  const Archive& operator<<(const std::optional<T>& v)
  {
    bool has_value = v.has_value();
    *this << has_value;
    if (has_value)
    {
      *this << v.value();
    }
    return *this;
  }

  template <class T, size_t size>
  Archive& operator>>(std::array<T, size>& v)
  {
    for (size_t i = 0; i < size; ++i)
    {
      T element;
      *this >> element;
      v[i] = std::move(element);
    }
    return *this;
  }

  template <typename T, size_t size>
  Archive& operator<<(const std::array<T, size>& v)
  {
    for (size_t i = 0; i < size; ++i)
    {
      *this << v[i];
    }
    return *this;
  }

  template <class T>
  Archive& operator>>(std::vector<T>& v)
  {
    size_t len;
    *this >> len;
    v.clear();
    v.reserve(len);
    for (size_t i = 0; i < len; ++i)
    {
      T element;
      *this >> element;
      v.push_back(std::move(element));
    }
    return *this;
  }

  template <class T>
  Archive& operator<<(const std::vector<T>& v)
  {
    size_t len = v.size();
    *this << len;
    for (const T& element : v)
    {
      *this << element;
    }
    return *this;
  }

  template <class T1, class T2>
  Archive& operator>>(std::map<T1, T2>& v)
  {
    size_t len;
    *this >> len;
    for (size_t i = 0; i < len; ++i)
    {
      std::pair<T1, T2> value;
      *this >> value;
      v[value.first] = value.second;
      // v.push_back(value);
    }
    return *this;
  }

  template <class T1, class T2>
  Archive& operator<<(const std::map<T1, T2>& v)
  {
    size_t len = v.size();
    *this << len;
    for (typename std::map<T1, T2>::const_iterator it = v.begin(); it != v.end(); ++it)
    {
      *this << *it;
    }
    return *this;
  }

  template <class T1, class T2>
  Archive& operator>>(std::pair<T1, T2>& v)
  {
    *this >> v.first >> v.second;
    return *this;
  }

  template <class T1, class T2>
  Archive& operator<<(const std::pair<T1, T2>& v)
  {
    *this << v.first << v.second;
    return *this;
  }

  Archive& operator>>(std::string& v)
  {
    size_t len;
    *this >> len;
    v.clear();
    char buffer[256];
    size_t toRead = len;
    while (toRead != 0)
    {
      size_t l = std::min(toRead, sizeof(buffer));
      stream.read(buffer, static_cast<std::streamsize>(l));
      if (!stream) throw std::runtime_error("malformed data");
      v += std::string(buffer, l);
      toRead -= l;
    }
    return *this;
  }

  Archive& operator<<(const std::string& v)
  {
    size_t len = v.length();
    *this << len;
    stream.write(v.c_str(), static_cast<std::streamsize>(len));
    return *this;
  }

 private:
  STREAM& stream;
};
