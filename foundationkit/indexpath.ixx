export module indexpath;

import <vector>;

export class IndexPath
{
public:
    IndexPath();
    IndexPath(const size_t index);
    size_t& operator[] (const size_t index);
    const size_t& operator[] (const size_t index) const;
    inline size_t lastIndex() const { if (!_path.empty()) return _path.back(); return 0; }
    const IndexPath operator+(const IndexPath& rhs);
    void increaseValueAtLastIndex();
    void decreaseValueAtLastIndex();
    size_t size();
    bool empty() const { return _path.empty(); }
    IndexPath appending(size_t index);
    IndexPath removingLastIndex() const;

    bool operator<(const IndexPath& otherObject) const;
    bool operator>(const IndexPath& otherObject) const;
    bool operator==(const IndexPath& otherObject) const;
private:
    std::vector<size_t> _path;
};
