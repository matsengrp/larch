#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

template <typename... Fs>
struct ElementStorage {
 public:
  using ExtraStorage = std::tuple<ExtraFeatureStorage<Fs>...>;

  template <typename Feature>
  static const bool contains_element_feature;

  template <typename CRTP>
  struct ConstElementViewBase : FeatureConstView<Fs, CRTP>... {};
  template <typename CRTP>
  struct MutableElementViewBase : ConstElementViewBase<CRTP>,
                                  FeatureMutableView<Fs, CRTP>... {
    using FeatureMutableView<Fs, CRTP>::operator=...;
  };

  ElementStorage() = default;
  MOVE_ONLY(ElementStorage);

  template <typename F>
  auto& GetFeatureStorage();

  template <typename F>
  const auto& GetFeatureStorage() const;

 private:
  std::tuple<Fs...> features_storage_;
};