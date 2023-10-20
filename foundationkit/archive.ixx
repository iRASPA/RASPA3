export module archive;

import <string>;
import <tuple>;
import <vector>;
import <map>;
import <istream>;
import <ostream>;
import <bit>;

export template <class STREAM> class Archive
{
public:
	Archive(STREAM &stream): stream(stream) {};

  template <class T> Archive& operator>>(T& v) \
  {
    stream.read(std::bit_cast<char*>(&v), sizeof(T)); 
    if(!stream) { throw std::runtime_error("malformed data"); } 
    if constexpr (std::endian::native == std::endian::little)
    {
      v = std::byteswap(v); 
    }
    return *this; 
  }

  template <class T> const Archive& operator<<(T v) const 
  {
    if constexpr (std::endian::native == std::endian::little)
    {
      v = std::byteswap(v);
    }
    stream.write(std::bit_cast<const char*>(&v), sizeof(T));
    return *this;
  }

  //a function that can serialize any enum
  template<typename Enum,
           typename = typename std::enable_if<std::is_enum<Enum>::value>::type>
  Archive& operator<<(const Enum& e) 
  {
    *this << static_cast<int64_t>(e);
    return stream;
  }

  //a function that can deserialize any enum from an archive
  template<typename Enum,
           typename = typename std::enable_if<std::is_enum<Enum>::value>::type>
  Archive& operator>>(Enum& e)
  {
    int64_t v;
    *this >> v;
    e = static_cast<Enum>(v);
    return stream;
  }


  template <class T> Archive& operator>>(std::vector<T>& v)
  {
    size_t len;
    *this >> len;
    for(size_t i = 0; i < len; ++i)
    {
      T element;
      *this >> element;
      v.push_back(element);
    }
    return *this;
  } 

  template <class T> const Archive& operator<<(const std::vector<T>& v) const
  { 
    size_t len = v.size();
    *this << len; 
    for(const T& element : v)
    {
      *this << element;
    }
    return *this;
  }

  template <class T1, class T2> Archive& operator>>(std::map<T1, T2>& v)
  {
    size_t len;
    *this >> len;
    for(size_t i = 0; i < len; ++i)
    { 
      std::pair<T1, T2> value;
      *this >> value;
      v.push_back(value);
    }
    return *this;
  }

  template <class T1, class T2> const Archive& operator<<(const std::map<T1, T2>& v) const
  {
    size_t len = v.size();
    *this << len;
    for(typename std::map<T1, T2>::const_iterator it = v.begin(); it != v.end(); ++it)
    {
      *this << *it;
    }
    return *this;
  }

  template <class T1, class T2> Archive& operator>>(std::pair<T1, T2>& v)
  {
    *this >> v.first >> v.second;
    return *this;
  }

  template <class T1, class T2> const Archive& operator<<(const std::pair<T1, T2>& v) const
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
    while(toRead != 0)
    {
      size_t l = std::min(toRead, sizeof(buffer));
      stream.read(buffer, l);
      if(!stream) throw std::runtime_error("malformed data");
      v += std::string(buffer, l);
      toRead -= l;
    }
    return *this;
  }

  const Archive& operator<<(const std::string& v) const
  {
    size_t len = v.length();
    *this << len;
    stream.write(v.c_str(), len);
    return *this;
  }

private:
  STREAM& stream;
};

