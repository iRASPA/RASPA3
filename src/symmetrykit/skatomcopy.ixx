module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <memory>
#include <vector>
#endif

export module skatom:skatomcopy;

#ifndef USE_LEGACY_HEADERS
import <memory>;
import <vector>;
#endif

import double3;

export class SKAtomCopy
{
 public:
  enum class AtomCopyType : int64_t
  {
    copy = 2,
    duplicate = 3
  };

  SKAtomCopy() : _position(), _tag(0), _type(AtomCopyType::copy) {}
  SKAtomCopy(const SKAtomCopy& atomCopy);

  double3 position() const { return _position; }
  void setPosition(double3 p) { _position = p; }
  AtomCopyType type() { return _type; }
  void setType(AtomCopyType type) { _type = type; }
  int64_t tag() { return _tag; }
  void setTag(int64_t tag) { _tag = tag; }
  int64_t asymmetricIndex() { return _asymmetricIndex; }
  void setAsymmetricIndex(int64_t value) { _asymmetricIndex = value; }

 private:
  [[maybe_unused]] int64_t _versionNumber{1};

  struct Hash
  {
    template <typename T>
    std::size_t operator()(T* const& p) const
    {
      return std::hash<T>()(*p);
    }
  };
  struct Compare
  {
    template <typename T>
    size_t operator()(T* const& a, T* const& b) const
    {
      return *a == *b;
    }
  };
  double3 _position;
  int64_t _tag;
  AtomCopyType _type;
  int64_t _asymmetricIndex;
};

SKAtomCopy::SKAtomCopy(const SKAtomCopy& atomCopy)
{
  this->_position = atomCopy._position;
  this->_type = atomCopy._type;
  this->_tag = atomCopy._tag;
  this->_asymmetricIndex = atomCopy._asymmetricIndex;
  // this->_bonds = {};
}
