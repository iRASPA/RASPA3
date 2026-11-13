module;

#if !defined(_WIN32)
#include <fcntl.h>
#include <unistd.h>
#endif

export module archive;

import std;

// on linux uint64_t is unsigned long        8
//          size_t   is unsigned long        8
//          size_t is an alias for uint64_t
// on mac   uint64_t is unsigned long long   8
//          size_t   is unsigned long        8
//          size_t is an alias for unsigned long

// The archive format uses the native endianness and native in-memory layout of the scalar types:
// restart files are scratch data for the machine that wrote them, and skipping per-element byte
// swapping makes checkpoint writes of large arrays a single bulk memcpy to the stream.

// Element types whose element-wise serialization is bit-identical to their in-memory layout, so
// contiguous containers of them can be read/written in one bulk stream operation. Enums are
// excluded (they are serialized widened to int64), bool is excluded (std::vector<bool> is packed
// and bool has no guaranteed representation), std::complex<double> is guaranteed to be laid out
// as {real, imag}.
// A type may declare that its archive form is its memory image ('static constexpr bool archiveIsMemoryImage'):
// contiguous containers of it are then streamed in bulk as well.
template <typename T>
concept ArchivedAsMemoryImage = requires {
  { T::archiveIsMemoryImage } -> std::convertible_to<bool>;
} && T::archiveIsMemoryImage;

template <typename T>
constexpr bool is_bulk_serializable_v = (std::is_arithmetic_v<T> && !std::is_same_v<T, bool>) ||
                                        std::is_same_v<T, std::complex<double>> || ArchivedAsMemoryImage<T>;

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

  // all arithmetic scalars (integers and floating point) are written in native byte order
  template <typename T>
    requires(std::is_arithmetic_v<T> && !std::is_same_v<T, bool>)
  Archive& operator>>(T& v)
  {
    stream.read(std::bit_cast<char*>(&v), sizeof(T));
    if (!stream)
    {
      throw std::runtime_error("malformed data");
    }
    return *this;
  }

  template <typename T>
    requires(std::is_arithmetic_v<T> && !std::is_same_v<T, bool>)
  Archive& operator<<(const T& v)
  {
    stream.write(std::bit_cast<const char*>(&v), sizeof(T));
    return *this;
  }

  Archive& operator>>(std::chrono::duration<double>& v)
  {
    double count{0.0};
    *this >> count;
    v = std::chrono::duration<double>(count);
    return *this;
  }

  Archive& operator<<(const std::chrono::duration<double>& v)
  {
    *this << v.count();
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

  Archive& operator<<(const std::complex<double>& v)
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
  Archive& operator<<(const std::optional<T>& v)
  {
    bool has_value = v.has_value();
    *this << has_value;
    if (has_value)
    {
      *this << v.value();
    }
    return *this;
  }

  template <class T1, class T2>
  Archive& operator>>(std::variant<T1, T2>& v)
  {
    std::size_t index{};
    *this >> index;
    if(index == 0)
    {
      T1 value{};
      *this >> value;
      v = std::move(value);
    }
    if(index == 1)
    {
      T2 value{};
      *this >> value;
      v = std::move(value);
    }
    // any other index is not a variant written by operator<<: the stream is out of step with the data
    if (index > 1) throw std::runtime_error("malformed data");
    return *this;
  }

  template <class T1, class T2>
  Archive& operator<<(const std::variant<T1, T2>& v)
  {
    std::size_t index = v.index();
    *this << index;
    if(index == 0)
    {
      *this << get<0>(v);
    }
    if(index == 1)
    {
      *this << get<1>(v);
    }
    return *this;
  }


  template <class T, std::size_t size>
  Archive& operator>>(std::array<T, size>& v)
  {
    if constexpr (is_bulk_serializable_v<T>)
    {
      stream.read(std::bit_cast<char*>(v.data()), static_cast<std::streamsize>(size * sizeof(T)));
      if (!stream)
      {
        throw std::runtime_error("malformed data");
      }
    }
    else
    {
      for (std::size_t i = 0; i < size; ++i)
      {
        T element;
        *this >> element;
        v[i] = std::move(element);
      }
    }
    return *this;
  }

  template <typename T, std::size_t size>
  Archive& operator<<(const std::array<T, size>& v)
  {
    if constexpr (is_bulk_serializable_v<T>)
    {
      stream.write(std::bit_cast<const char*>(v.data()), static_cast<std::streamsize>(size * sizeof(T)));
    }
    else
    {
      for (std::size_t i = 0; i < size; ++i)
      {
        *this << v[i];
      }
    }
    return *this;
  }

  template <class T>
  Archive& operator>>(std::vector<T>& v)
  {
    std::size_t len;
    *this >> len;
    if constexpr (is_bulk_serializable_v<T>)
    {
      v.resize(len);
      stream.read(std::bit_cast<char*>(v.data()), static_cast<std::streamsize>(len * sizeof(T)));
      if (!stream)
      {
        throw std::runtime_error("malformed data");
      }
    }
    else
    {
      v.clear();
      v.reserve(len);
      for (std::size_t i = 0; i < len; ++i)
      {
        T element;
        *this >> element;
        v.push_back(std::move(element));
      }
    }
    return *this;
  }

  template <class T>
  Archive& operator<<(const std::vector<T>& v)
  {
    std::size_t len = v.size();
    *this << len;
    if constexpr (is_bulk_serializable_v<T>)
    {
      stream.write(std::bit_cast<const char*>(v.data()), static_cast<std::streamsize>(len * sizeof(T)));
    }
    else
    {
      for (const T& element : v)
      {
        *this << element;
      }
    }
    return *this;
  }

  template <class T1, class T2>
  Archive& operator>>(std::map<T1, T2>& v)
  {
    std::size_t len;
    *this >> len;
    // the archive holds the whole map: what was in it before is replaced, not merged
    v.clear();
    for (std::size_t i = 0; i < len; ++i)
    {
      std::pair<T1, T2> value;
      *this >> value;
      v[std::move(value.first)] = std::move(value.second);
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
  Archive& operator>>(std::unordered_map<T1, T2>& v)
  {
    std::size_t len;
    *this >> len;
    // the archive holds the whole map: what was in it before is replaced, not merged
    v.clear();
    v.reserve(len);
    for (std::size_t i = 0; i < len; ++i)
    {
      std::pair<T1, T2> value;
      *this >> value;
      v[std::move(value.first)] = std::move(value.second);
    }
    return *this;
  }

  template <class T1, class T2>
  Archive& operator<<(const std::unordered_map<T1, T2>& v)
  {
    std::size_t len = v.size();
    *this << len;
    for (typename std::unordered_map<T1, T2>::const_iterator it = v.begin(); it != v.end(); ++it)
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
    // sized once and read in one call (it was appended in 256-byte pieces through a temporary string each)
    v.resize(len);
    if (len != 0) stream.read(v.data(), static_cast<std::streamsize>(len));
    if (!stream) throw std::runtime_error("malformed data");
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

// Byte-order marker at the start of every binary restart file. The archive format uses the
// native byte order, so a file written on a machine with different endianness reads this back
// permuted and can be rejected with a clear error instead of deserializing garbage.
constexpr std::uint32_t restartFileByteOrderMarker = 0x01020304;

// Restart file layout: byte-order marker, payload size, payload checksum, payload.
constexpr std::uint64_t restartFileHeaderSize = sizeof(std::uint32_t) + 2uz * sizeof(std::uint64_t);

// Streaming 64-bit hash over the payload (8-byte words with splitmix64-style mixing, seeded with
// the payload size). Not cryptographic; detects truncation and bit corruption. Throws when the
// stream ends before 'payloadSize' bytes.
// The checksum of a restart payload, taken over the bytes as a stream: eight-byte words in file order (the last
// one zero-padded), each folded into the running value, and the payload size folded in last. Fed by the writer
// as the payload goes out (through 'HashingStreambuf') and by the reader as it comes in, so a file is written
// once and never read back to be summed. (The size used to seed the hash, which made a single pass over an
// outgoing payload impossible; a checkpoint from before this change fails its checksum and is rejected, as a
// corrupt one is.)
struct RestartPayloadHash
{
  static constexpr std::uint64_t mixBits(std::uint64_t h)
  {
    h ^= h >> 30;
    h *= 0xbf58476d1ce4e5b9;
    h ^= h >> 27;
    h *= 0x94d049bb133111eb;
    h ^= h >> 31;
    return h;
  }

  std::uint64_t hash{0xcbf29ce484222325};
  std::uint64_t size{0};
  std::array<char, 8> pending{};  // the bytes of a word not yet complete across calls
  std::size_t pendingBytes{0};

  void update(const char* data, std::size_t count)
  {
    size += count;
    // complete a word left pending by the previous call
    while (count != 0 && pendingBytes != 0)
    {
      pending[pendingBytes++] = *data++;
      --count;
      if (pendingBytes == 8)
      {
        std::uint64_t word{0};
        std::memcpy(&word, pending.data(), 8);
        hash = mixBits(hash ^ word);
        pendingBytes = 0;
      }
    }
    for (; count >= 8; data += 8, count -= 8)
    {
      std::uint64_t word{0};
      std::memcpy(&word, data, 8);
      hash = mixBits(hash ^ word);
    }
    for (; count != 0; --count) pending[pendingBytes++] = *data++;
  }

  [[nodiscard]] std::uint64_t finish() const
  {
    std::uint64_t h = hash;
    if (pendingBytes != 0)
    {
      std::uint64_t word{0};
      std::memcpy(&word, pending.data(), pendingBytes);
      h = mixBits(h ^ word);
    }
    return mixBits(mixBits(h) ^ size);
  }
};

// A stream buffer that hands every byte to a RestartPayloadHash on its way to the underlying buffer.
class HashingStreambuf : public std::streambuf
{
 public:
  HashingStreambuf(std::streambuf* sink, RestartPayloadHash& hash) : sink_(sink), hash_(hash) {}

 protected:
  std::streamsize xsputn(const char* data, std::streamsize count) override
  {
    const std::streamsize written = sink_->sputn(data, count);
    if (written > 0) hash_.update(data, static_cast<std::size_t>(written));
    return written;
  }
  int overflow(int character) override
  {
    if (character == traits_type::eof())
    {
      return sink_->pubsync() == 0 ? traits_type::not_eof(character) : traits_type::eof();
    }
    const char byte = traits_type::to_char_type(character);
    if (sink_->sputc(byte) == traits_type::eof()) return traits_type::eof();
    hash_.update(&byte, 1);
    return character;
  }
  int sync() override { return sink_->pubsync(); }

 private:
  std::streambuf* sink_;
  RestartPayloadHash& hash_;
};

// The checksum of the 'payloadSize' bytes that follow in 'stream' (the reader's side). Throws when the stream
// ends before 'payloadSize' bytes.
inline std::uint64_t hashRestartFilePayload(std::istream& stream, std::uint64_t payloadSize)
{
  RestartPayloadHash hash{};
  std::vector<char> buffer(1uz << 20uz);
  std::uint64_t remaining = payloadSize;
  while (remaining != 0uz)
  {
    const std::size_t chunk = static_cast<std::size_t>(std::min<std::uint64_t>(remaining, buffer.size()));
    stream.read(buffer.data(), static_cast<std::streamsize>(chunk));
    if (!stream)
    {
      throw std::runtime_error("truncated restart file payload");
    }
    hash.update(buffer.data(), chunk);
    remaining -= chunk;
  }
  return hash.finish();
}

// Durability: flush the written file contents to stable storage before the rename, so a power
// loss cannot leave a fully-renamed but partially-persisted restart file. On macOS fsync only
// flushes to the drive cache; F_FULLFSYNC is the real storage barrier. Returns false when the data
// could not be persisted (a deferred write error such as ENOSPC or EIO surfaces here), so the file
// must not be published.
inline bool syncFileToStableStorage(const std::string& fileName) noexcept
{
#if !defined(_WIN32)
  int fileDescriptor = ::open(fileName.c_str(), O_RDONLY);
  if (fileDescriptor < 0) return false;
#if defined(__APPLE__) && defined(__MACH__)
  bool synced = ::fcntl(fileDescriptor, F_FULLFSYNC) == 0 || ::fsync(fileDescriptor) == 0;
#else
  bool synced = ::fsync(fileDescriptor) == 0;
#endif
  if (::close(fileDescriptor) != 0) synced = false;
  return synced;
#else
  return true;
#endif
}

// Durability: persist the renames themselves (directory entries) to stable storage.
inline void syncParentDirectory(const std::string& fileName) noexcept
{
#if !defined(_WIN32)
  std::filesystem::path parent = std::filesystem::path(fileName).parent_path();
  if (parent.empty()) parent = ".";
  int fileDescriptor = ::open(parent.c_str(), O_RDONLY);
  if (fileDescriptor < 0) return;
  ::fsync(fileDescriptor);
  ::close(fileDescriptor);
#endif
}

/**
 * \brief Serializes 'object' to a binary restart file, atomically replacing any previous file.
 *
 * The object is written to '<fileName>_temp' using a large stream buffer (few, large write
 * system calls), the header is completed with the payload size and checksum, and the file is
 * flushed to stable storage before the renames, so neither a crash nor a power loss can destroy
 * the previous valid restart file. The previous file is kept as '<fileName>.prev' so one older
 * checkpoint always remains as a fallback against corruption of the latest file. Failures leave
 * the previous restart file intact and are not propagated: a failed periodic checkpoint should
 * not kill a running simulation.
 */
export template <typename T>
void writeBinaryRestartFile(const T& object, const std::string& fileName = "restart_data.bin")
{
  try
  {
    const std::string temporaryFileName = fileName + "_temp";

    // one pass: the header with placeholders, the payload summed as it is written, the header patched in place
    {
      std::vector<char> buffer(1uz << 20uz);
      std::ofstream ofile;
      ofile.rdbuf()->pubsetbuf(buffer.data(), static_cast<std::streamsize>(buffer.size()));
      ofile.open(temporaryFileName, std::ios::binary);

      Archive<std::ofstream> archive(ofile);
      archive << restartFileByteOrderMarker;
      archive << std::uint64_t{0};  // payload size, patched below
      archive << std::uint64_t{0};  // payload checksum, patched below

      // the payload goes out through a buffer that sums it on the way
      std::filebuf* fileBuffer = ofile.rdbuf();
      RestartPayloadHash hash{};
      HashingStreambuf hashingBuffer(fileBuffer, hash);
      ofile.std::ios::rdbuf(&hashingBuffer);
      archive << object;
      ofile.flush();
      // rdbuf(sb) clears the stream state: carry a failed payload write (ENOSPC, quota, EIO) across the swap, or
      // the header below would seal the truncated payload with a valid checksum
      const std::ios::iostate payloadState = ofile.rdstate();
      ofile.std::ios::rdbuf(fileBuffer);
      ofile.setstate(payloadState);
      if (!ofile)
      {
        ofile.close();
        std::error_code error;
        std::filesystem::remove(temporaryFileName, error);
        return;
      }

      const std::uint64_t payloadSize = hash.size;
      const std::uint64_t checksum = hash.finish();
      ofile.seekp(static_cast<std::streamoff>(sizeof(std::uint32_t)));
      ofile.write(std::bit_cast<const char*>(&payloadSize), sizeof(std::uint64_t));
      ofile.write(std::bit_cast<const char*>(&checksum), sizeof(std::uint64_t));
      ofile.close();
      if (!ofile)
      {
        std::error_code error;
        std::filesystem::remove(temporaryFileName, error);
        return;
      }
    }

    if (!syncFileToStableStorage(temporaryFileName))
    {
      std::error_code error;
      std::filesystem::remove(temporaryFileName, error);
      return;
    }

    // keep the previous checkpoint as fallback, then atomically publish the new one
    std::error_code error;
    if (std::filesystem::exists(fileName, error))
    {
      std::filesystem::rename(fileName, fileName + ".prev", error);
    }
    std::filesystem::rename(temporaryFileName, fileName, error);

    syncParentDirectory(fileName);
  }
  catch (const std::exception&)
  {
  }
}

// Reads and fully validates a single restart file: byte-order marker, payload size against the
// actual file size (truncation), and payload checksum (bit corruption) before deserializing.
template <typename T>
void readValidatedBinaryRestartFile(T& object, const std::string& fileName)
{
  std::ifstream ifile(fileName, std::ios::binary);
  if (!ifile.is_open())
  {
    throw std::runtime_error(std::format("Binary restart file '{}' doesn't exist\n", fileName));
  }

  // a file shorter than the header cannot hold the sizes read below (and the payload size would underflow)
  const std::uint64_t fileSize = static_cast<std::uint64_t>(std::filesystem::file_size(fileName));
  if (fileSize < restartFileHeaderSize)
  {
    throw std::runtime_error(
        std::format("Binary restart file '{}' is truncated ({} bytes, shorter than the "
                    "{}-byte header)\n",
                    fileName, fileSize, restartFileHeaderSize));
  }

  Archive<std::ifstream> archive(ifile);
  std::uint32_t byteOrderMarker;
  archive >> byteOrderMarker;
  if (byteOrderMarker != restartFileByteOrderMarker)
  {
    throw std::runtime_error(
        std::format("Binary restart file '{}' was written on a machine with an incompatible "
                    "byte order (or is not a restart file)\n",
                    fileName));
  }

  std::uint64_t storedPayloadSize, storedChecksum;
  archive >> storedPayloadSize;
  archive >> storedChecksum;

  const std::uint64_t actualPayloadSize = fileSize - restartFileHeaderSize;
  if (actualPayloadSize != storedPayloadSize)
  {
    throw std::runtime_error(std::format("Binary restart file '{}' is truncated or padded ({} payload "
                                         "bytes, expected {})\n",
                                         fileName, actualPayloadSize, storedPayloadSize));
  }

  const std::uint64_t checksum = hashRestartFilePayload(ifile, storedPayloadSize);
  if (checksum != storedChecksum)
  {
    throw std::runtime_error(std::format("Binary restart file '{}' is corrupt (checksum mismatch)\n", fileName));
  }

  ifile.clear();
  ifile.seekg(static_cast<std::streamoff>(restartFileHeaderSize));
  archive >> object;
}

/**
 * \brief Deserializes 'object' from a binary restart file written by writeBinaryRestartFile.
 *
 * The file is fully validated (byte order, size, checksum) before deserialization. When the
 * latest file is missing or corrupt and a previous checkpoint '<fileName>.prev' exists, that
 * fallback is used instead (with a warning): resuming from one checkpoint earlier beats losing
 * the whole run.
 */
export template <typename T>
void readBinaryRestartFile(T& object, const std::string& fileName = "restart_data.bin")
{
  const std::string previousFileName = fileName + ".prev";

  try
  {
    readValidatedBinaryRestartFile(object, fileName);
    return;
  }
  catch (const std::exception& error)
  {
    if (!std::filesystem::exists(previousFileName))
    {
      throw;
    }
    std::cerr << std::format("Warning: {}Falling back to the previous checkpoint '{}'\n", error.what(),
                             previousFileName);
  }

  readValidatedBinaryRestartFile(object, previousFileName);
}
