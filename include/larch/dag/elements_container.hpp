#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

/**
 * Stores a collection of elements (nodes or edges, distinguished by the `Id`
 * parameter).
 */
template <typename Id, typename ElementStorageT, typename... Features>
struct ElementsContainer {
 public:
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

  size_t GetCount() const;

  Id Append();

  void Add(Id id);

  void Initialize(size_t size);

  template <typename Feature>
  auto& GetFeatureStorage(Id id);
  template <typename Feature>
  const auto& GetFeatureStorage(Id id) const;

  template <typename Feature>
  auto& GetFeatureExtraStorage();
  template <typename Feature>
  const auto& GetFeatureExtraStorage() const;

 private:
  std::vector<ElementStorageT> elements_storage_;
  std::vector<std::tuple<Features...>> features_storage_;
  std::tuple<ExtraFeatureStorage<Features>...> extra_features_storage_;
  typename ElementStorageT::ExtraStorage elements_extra_features_storage_;
};