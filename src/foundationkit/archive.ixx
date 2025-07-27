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
import std;
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
    *this << static_cast<std::int64_t>(e);
    return *this;
  }

  // a function that can deserialize any enum from an archive
  template <typename Enum, typename = typename std::enable_if<std::is_enum<Enum>::value>::type>
  Archive& operator>>(Enum& e)
  {
    std::int64_t v;
    *this >> v;
    e = static_cast<Enum>(v);
    return *this;
  }
  Archive& operator>>(bool& v)
  {
    std::uint8_t w;
    *this >> w;
    v = static_cast<bool>(w);
    return *this;
  }

  Archive& operator<<(const bool& v)
  {
    *this << static_cast<std::uint8_t>(v);
    return *this;
  }

  Archive& operator>>(std::int8_t& v)
  {
    stream.read(std::bit_cast<char*>(&v), sizeof(std::int8_t));
    if (!stream)
    {
      throw std::runtime_error("malformed data");
    }
    return *this;
  }

  Archive& operator<<(const std::int8_t& v)
  {
    stream.write(std::bit_cast<const char*>(&v), sizeof(std::int8_t));
    return *this;
  }

  Archive& operator>>(std::uint8_t& v)
  {
    stream.read(std::bit_cast<char*>(&v), sizeof(std::uint8_t));
    if (!stream)
    {
      throw std::runtime_error("malformed data");
    }
    return *this;
  }

  Archive& operator<<(const std::uint8_t& v)
  {
    stream.write(std::bit_cast<const char*>(&v), sizeof(std::uint8_t));
    return *this;
  }

  Archive& operator>>(std::int16_t& v)
  {
    stream.read(std::bit_cast<char*>(&v), sizeof(std::int16_t));
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

  Archive& operator<<(const std::int16_t& v)
  {
    std::int16_t w{v};
    if constexpr (std::endian::native == std::endian::little)
    {
      w = std::byteswap(w);
    }
    stream.write(std::bit_cast<const char*>(&w), sizeof(std::int16_t));
    return *this;
  }

  Archive& operator>>(std::uint16_t& v)
  {
    stream.read(std::bit_cast<char*>(&v), sizeof(std::uint16_t));
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

  Archive& operator<<(const std::uint16_t& v)
  {
    std::uint16_t w{v};
    if constexpr (std::endian::native == std::endian::little)
    {
      w = std::byteswap(w);
    }
    stream.write(std::bit_cast<const char*>(&w), sizeof(std::uint16_t));
    return *this;
  }

  Archive& operator>>(std::int32_t& v)
  {
    stream.read(std::bit_cast<char*>(&v), sizeof(std::int32_t));
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

  Archive& operator<<(const std::int32_t& v)
  {
    std::int32_t w{v};
    if constexpr (std::endian::native == std::endian::little)
    {
      w = std::byteswap(w);
    }
    stream.write(std::bit_cast<const char*>(&w), sizeof(std::int32_t));
    return *this;
  }

  Archive& operator>>(std::uint32_t& v)
  {
    stream.read(std::bit_cast<char*>(&v), sizeof(std::uint32_t));
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

  Archive& operator<<(const std::uint32_t& v)
  {
    std::uint32_t w{v};
    if constexpr (std::endian::native == std::endian::little)
    {
      w = std::byteswap(w);
    }
    stream.write(std::bit_cast<const char*>(&w), sizeof(std::uint32_t));
    return *this;
  }

  Archive& operator>>(std::int64_t& v)
  {
    stream.read(std::bit_cast<char*>(&v), sizeof(std::int64_t));
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

  Archive& operator<<(const std::int64_t& v)
  {
    std::int64_t w{v};
    if constexpr (std::endian::native == std::endian::little)
    {
      w = std::byteswap(w);
    }
    stream.write(std::bit_cast<const char*>(&w), sizeof(std::int64_t));
    return *this;
  }

#if defined(__APPLE__) && defined(__MACH__)
  Archive& operator>>(std::uint64_t& v)
  {
    stream.read(std::bit_cast<char*>(&v), sizeof(std::uint64_t));
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

  Archive& operator<<(const std::uint64_t& v)
  {
    std::uint64_t w{v};
    if constexpr (std::endian::native == std::endian::little)
    {
      w = std::byteswap(w);
    }
    stream.write(std::bit_cast<const char*>(&w), sizeof(std::uint64_t));
    return *this;
  }
#endif

  Archive& operator>>(std::size_t& v)
  {
    stream.read(std::bit_cast<char*>(&v), sizeof(std::size_t));
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

  Archive& operator<<(const std::size_t& v)
  {
    std::size_t w{v};
    if constexpr (std::endian::native == std::endian::little)
    {
      w = std::byteswap(w);
    }
    stream.write(std::bit_cast<const char*>(&w), sizeof(std::size_t));
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

  template <class T, std::size_t size>
  Archive& operator>>(std::array<T, size>& v)
  {
    for (std::size_t i = 0; i < size; ++i)
    {
      T element;
      *this >> element;
      v[i] = std::move(element);
    }
    return *this;
  }

  template <typename T, std::size_t size>
  Archive& operator<<(const std::array<T, size>& v)
  {
    for (std::size_t i = 0; i < size; ++i)
    {
      *this << v[i];
    }
    return *this;
  }

  template <class T>
  Archive& operator>>(std::vector<T>& v)
  {
    std::size_t len;
    *this >> len;
    v.clear();
    v.reserve(len);
    for (std::size_t i = 0; i < len; ++i)
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
    std::size_t len = v.size();
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
    std::size_t len;
    *this >> len;
    for (std::size_t i = 0; i < len; ++i)
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
    std::size_t len = v.size();
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

  template <class T1, class T2, class T3>
  Archive& operator>>(std::tuple<T1, T2, T3>& v)
  {
    *this >> std::get<0>(v) >> std::get<1>(v) >> std::get<2>(v);
    return *this;
  }

  template <class T1, class T2, class T3>
  Archive& operator<<(const std::tuple<T1, T2, T3>& v)
  {
    *this << std::get<0>(v) << std::get<1>(v) << std::get<2>(v);
    return *this;
  }

  template <class T1, class T2, class T3, class T4>
  Archive& operator>>(std::tuple<T1, T2, T3, T4>& v)
  {
    *this >> std::get<0>(v) >> std::get<1>(v) >> std::get<2>(v) >> std::get<3>(v);
    return *this;
  }

  template <class T1, class T2, class T3, class T4>
  Archive& operator<<(const std::tuple<T1, T2, T3, T4>& v)
  {
    *this << std::get<0>(v) << std::get<1>(v) << std::get<2>(v) << std::get<3>(v);
    return *this;
  }

  Archive& operator>>(std::string& v)
  {
    std::size_t len;
    *this >> len;
    v.clear();
    char buffer[256];
    std::size_t toRead = len;
    while (toRead != 0)
    {
      std::size_t l = std::min(toRead, sizeof(buffer));
      stream.read(buffer, static_cast<std::streamsize>(l));
      if (!stream) throw std::runtime_error("malformed data");
      v += std::string(buffer, l);
      toRead -= l;
    }
    return *this;
  }

  Archive& operator<<(const std::string& v)
  {
    std::size_t len = v.length();
    *this << len;
    stream.write(v.c_str(), static_cast<std::streamsize>(len));
    return *this;
  }

 private:
  STREAM& stream;
};
