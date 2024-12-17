#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

/**
 * Stores a collection of elements (nodes or edges, distinguished by the `Id`
 * parameter).
 */
template <Component C, typename ElementStorageT,
          IdContinuity IdCont = IdContinuity::Dense, typename... Features>
struct ElementsContainer {
 public:
  using FeatureTypes = std::tuple<Features...>;
  using AllFeatureTypes = decltype(std::tuple_cat(
      FeatureTypes{}, typename ElementStorageT::FeatureTypes{}));

  static constexpr IdContinuity id_continuity = IdCont;

  template <typename Feature>
  static const bool contains_element_feature;

  template <typename CRTP>
  struct ConstElementViewBase : ElementStorageT::template ConstElementViewBase<CRTP>,
                                FeatureConstView<Features, CRTP>... {};
  template <typename CRTP>
  struct MutableElementViewBase
      : ElementStorageT::template MutableElementViewBase<CRTP>,
        FeatureMutableView<Features, CRTP>... {
    using ElementStorageT::template MutableElementViewBase<CRTP>::operator=;
    using FeatureMutableView<Features, CRTP>::operator=...;
  };

  template <typename CRTP>
  struct ExtraConstElementViewBase
      : ElementStorageT::template ExtraConstElementViewBase<CRTP>,
        ExtraFeatureConstView<Features, CRTP>... {};
  template <typename CRTP>
  struct ExtraMutableElementViewBase
      : ElementStorageT::template ExtraMutableElementViewBase<CRTP>,
        ExtraFeatureMutableView<Features, CRTP>... {};

  ElementsContainer() = default;
  MOVE_ONLY(ElementsContainer);

  template <typename VT>
  size_t GetCount() const;

  template <typename VT>
  Id<C> GetNextAvailableId() const {
    return {GetCount<VT>()};
  }

  template <typename VT>
  bool ContainsId(Id<C> id) const {
    return id.value < GetCount<VT>();
  }

  Id<C> Append();

  void Add(Id<C> id);

  void Initialize(size_t size);

  void Clear();

  template <typename Feature>
  auto& GetFeatureStorage(Id<C> id);
  template <typename Feature>
  const auto& GetFeatureStorage(Id<C> id) const;

  template <typename Feature>
  auto& GetFeatureExtraStorage();
  template <typename Feature>
  const auto& GetFeatureExtraStorage() const;

  template <typename VT>
  auto All() const {
    return ranges::views::iota(size_t{0}, GetCount<VT>()) |
           ranges::views::transform([](size_t i) -> Id<C> { return {i}; });
  }

 private:
  IdContainer<Id<C>, ElementStorageT, id_continuity> elements_storage_;
  IdContainer<Id<C>, std::tuple<Features...>, id_continuity> features_storage_;
  std::tuple<ExtraFeatureStorage<Features>...> extra_features_storage_;
  typename ElementStorageT::ExtraStorage elements_extra_features_storage_;
};
