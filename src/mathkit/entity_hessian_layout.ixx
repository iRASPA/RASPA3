module;

export module entity_hessian_layout;

import std;

/**
 * Row-major layout for an asymmetric dense Hessian matrix.
 */
export struct layout_entity_hessian
{
  template <class Extents>
    requires(Extents::rank() == 2)
  struct mapping
  {
    using extents_type = Extents;
    using size_type = typename extents_type::size_type;

    std::size_t nDofs{};

    constexpr mapping() noexcept = default;

    constexpr mapping(const extents_type &extents) noexcept : nDofs(extents.extent(0))
    {
      if constexpr (extents_type::rank_dynamic() == 0)
      {
        static_assert(extents_type::static_extent(0) == extents_type::static_extent(1));
      }
    }

    constexpr const extents_type &extents() const noexcept { return _extents; }

    constexpr size_type required_span_size() const noexcept { return nDofs * nDofs; }

    constexpr size_type operator()(size_type i, size_type j) const noexcept { return i * nDofs + j; }

   private:
    extents_type _extents{};
  };
};
