module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <vector>
#endif

export module indexpath;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

export class IndexPath
{
 public:
  IndexPath();
  IndexPath(const std::size_t index);
  std::size_t& operator[](const std::size_t index);
  const std::size_t& operator[](const std::size_t index) const;
  inline std::size_t lastIndex() const
  {
    if (!_path.empty()) return _path.back();
    return 0;
  }
  const IndexPath operator+(const IndexPath& rhs);
  void increaseValueAtLastIndex();
  void decreaseValueAtLastIndex();
  std::size_t size();
  bool empty() const { return _path.empty(); }
  IndexPath appending(std::size_t index);
  IndexPath removingLastIndex() const;

  bool operator<(const IndexPath& otherObject) const;
  bool operator>(const IndexPath& otherObject) const;
  bool operator==(const IndexPath& otherObject) const;

 private:
  std::vector<std::size_t> _path;
};
