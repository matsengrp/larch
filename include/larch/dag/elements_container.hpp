#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

/**
 * Stores a collection of elements (nodes or edges, distinguished by the `Id`
 * parameter).
 */
template <typename Id, typename ES, typename... Fs>
struct ElementsContainer {
 public:
  template <typename Feature>
  static const bool contains_element_feature;

  template <typename CRTP>
  struct ConstElementViewBase : ES::template ConstElementViewBase<CRTP>,
                                FeatureConstView<Fs, CRTP>... {};
  template <typename CRTP>
  struct MutableElementViewBase : ES::template MutableElementViewBase<CRTP>,
                                  FeatureMutableView<Fs, CRTP>... {
    using ES::template MutableElementViewBase<CRTP>::operator=;
    using FeatureMutableView<Fs, CRTP>::operator=...;
  };

  ElementsContainer() = default;
  MOVE_ONLY(ElementsContainer);

  size_t GetCount() const;

  Id Append();

  void Add(Id id);

  void Initialize(size_t size);

  template <typename F>
  auto& GetFeatureStorage(Id id);
  template <typename F>
  const auto& GetFeatureStorage(Id id) const;

  template <typename F>
  auto& GetFeatureExtraStorage();
  template <typename F>
  const auto& GetFeatureExtraStorage() const;

 private:
  template <typename, typename, typename...>
  friend struct DAGStorage;  // TODO remove

  std::vector<ES> elements_storage_;
  std::vector<std::tuple<Fs...>> features_storage_;
  std::tuple<ExtraFeatureStorage<Fs>...> extra_features_storage_;
  typename ES::ExtraStorage elements_extra_features_storage_;
};