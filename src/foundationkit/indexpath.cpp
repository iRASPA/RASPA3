module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cstddef>
#include <vector>
#endif

module indexpath;

#ifndef USE_LEGACY_HEADERS
import <algorithm>;
import <vector>;
#endif

IndexPath::IndexPath() { _path.reserve(10); }

IndexPath::IndexPath(const size_t index)
{
  _path.reserve(10);
  _path.push_back(index);
}

size_t& IndexPath::operator[](const size_t index)  // for non-const objects: can be used for assignment
{
  return _path[index];
}

const size_t& IndexPath::operator[](const size_t index) const  // for const objects: can only be used for access
{
  return _path[index];
}

void IndexPath::increaseValueAtLastIndex()
{
  if (!_path.empty())
  {
    _path.back() += 1;
  }
}

void IndexPath::decreaseValueAtLastIndex()
{
  if (!_path.empty())
  {
    _path.back() -= 1;
  }
}

size_t IndexPath::size() { return _path.size(); }

IndexPath IndexPath::appending(size_t index)
{
  IndexPath indexpath{};
  indexpath._path.insert(indexpath._path.end(), _path.begin(), _path.end());
  indexpath._path.push_back(index);
  return indexpath;
}

IndexPath IndexPath::removingLastIndex() const
{
  IndexPath indexPath{};
  if (_path.size() <= 1)
  {
    return indexPath;
  }
  indexPath._path = std::vector<size_t>(_path.begin(), std::prev(_path.end(), 1));
  return indexPath;
}

const IndexPath IndexPath::operator+(const IndexPath& rhs)
{
  IndexPath indexpath{};
  indexpath._path.insert(indexpath._path.end(), _path.begin(), _path.end());
  indexpath._path.insert(indexpath._path.end(), rhs._path.begin(), rhs._path.end());
  return indexpath;
}

bool IndexPath::operator<(const IndexPath& otherObject) const
{
  size_t l1 = _path.size();
  size_t l2 = otherObject._path.size();
  for (size_t pos = 0; pos < std::min(l1, l2); pos++)
  {
    size_t i1 = _path[pos];
    size_t i2 = otherObject._path[pos];
    if (i1 < i2)
    {
      return true;
    }
    else if (i1 > i2)
    {
      return false;
    }
  }
  if (l1 < l2)
  {
    return true;
  }
  else if (l1 > l2)
  {
    return false;
  }
  return false;
}

bool IndexPath::operator>(const IndexPath& otherObject) const
{
  size_t l1 = _path.size();
  size_t l2 = otherObject._path.size();
  for (size_t pos = 0; pos < std::min(l1, l2); pos++)
  {
    size_t i1 = _path[pos];
    size_t i2 = otherObject._path[pos];
    if (i1 < i2)
    {
      return false;
    }
    else if (i1 > i2)
    {
      return true;
    }
  }
  if (l1 < l2)
  {
    return false;
  }
  else if (l1 > l2)
  {
    return true;
  }
  return false;
}

bool IndexPath::operator==(const IndexPath& otherObject) const
{
  size_t l1 = _path.size();
  size_t l2 = otherObject._path.size();
  for (size_t pos = 0; pos < std::min(l1, l2); pos++)
  {
    size_t i1 = _path[pos];
    size_t i2 = otherObject._path[pos];
    if (i1 < i2)
    {
      return false;
    }
    else if (i1 > i2)
    {
      return false;
    }
  }
  if (l1 < l2)
  {
    return false;
  }
  else if (l1 > l2)
  {
    return false;
  }
  return false;
}
