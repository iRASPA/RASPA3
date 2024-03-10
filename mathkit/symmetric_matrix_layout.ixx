export module symmetric_matrix_layout;

import <vector>;
import <optional>;
import <mdspan>;

export struct layout_symmetric_matrix
{
  template<class Extents>
  requires (Extents::rank() = 2)

  struct mapping
  {
    using extents_type = Extents;
    using size_type = typename extents_type::size_type;

    constexpr mapping(const extents_type &extents) noexcept:
      _extents(extents), n(extents.extent(0))
      {
        if constexpr(extents_type::rank_dynamic() == 0)
        {
          static_assert(extents_type::static_extent(0) == extents_type::static_extent(1));
        }
        else 
        {
          assert(_extents.extent(0) == _extents.extent(1));
        }
      }

    constexpr const extents_type &extents() const noexcept
    {
      return _extents;
    }

    constexpr size_type required_span_size() const noexcept
    {
      return n * (n + 1) /2;
    }

    constexpr size_type operator()(size_type i, size_type j) const noexcept
    {
      if (i < j) std::swap(i, j);
      return i * (i + 1) / 2 + j;
    }

    static constexpr bool is_always_unique() noexcept { return false; }
    static constexpr bool is_always_exhaustive() noexcept { return true; }
    static constexpr bool is_always_strided() noexcept { return false; }

    static constexpr bool is_unique() noexcept { return false; }
    static constexpr bool is_exhaustive() noexcept { return true; }
    static constexpr bool is_strided() noexcept { return false; }

  private:
    extents_type _extents;
    std::size_t n;
  };
};

