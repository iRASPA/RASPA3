module;

module indexpath;

import std;

IndexPath::IndexPath() {}

IndexPath::IndexPath(const std::size_t index)
{
  _path.reserve(10);
  _path.push_back(index);
}

std::size_t& IndexPath::operator[](const std::size_t index)  // for non-const objects: can be used for assignment
{
  return _path[index];
}

const std::size_t& IndexPath::operator[](
    const std::size_t index) const  // for const objects: can only be used for access
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
  // an index of zero stays zero: the subtraction would wrap to the largest std::size_t
  if (!_path.empty() && _path.back() > 0)
  {
    _path.back() -= 1;
  }
}

std::size_t IndexPath::size() const { return _path.size(); }

IndexPath IndexPath::appending(std::size_t index) const
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
  indexPath._path = std::vector<std::size_t>(_path.begin(), std::prev(_path.end(), 1));
  return indexPath;
}

const IndexPath IndexPath::operator+(const IndexPath& rhs) const
{
  IndexPath indexpath{};
  indexpath._path.insert(indexpath._path.end(), _path.begin(), _path.end());
  indexpath._path.insert(indexpath._path.end(), rhs._path.begin(), rhs._path.end());
  return indexpath;
}

bool IndexPath::operator<(const IndexPath& otherObject) const
{
  return std::lexicographical_compare(_path.begin(), _path.end(), otherObject._path.begin(), otherObject._path.end());
}

bool IndexPath::operator>(const IndexPath& otherObject) const
{
  return std::lexicographical_compare(otherObject._path.begin(), otherObject._path.end(), _path.begin(), _path.end());
}

bool IndexPath::operator==(const IndexPath& otherObject) const
{
  return _path == otherObject._path;
}
