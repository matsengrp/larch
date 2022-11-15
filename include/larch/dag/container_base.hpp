#include <type_traits>
#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

template <typename Id, template <typename> typename ViewType, typename Storage,
          typename... Features>
class DefaultContainerBase {
 public:
  using StorageType = Storage;
  using FeaturesType = std::tuple<std::vector<Features>...>;

  using FeaturesGlobalData =
      filter_t<std::is_class, std::tuple,
               typename FeatureTraits<
                   Features, ViewType<DAGView<Storage>>>::template GlobalData<Id>...>;

  template <typename Feature>
  using GlobalDataForFeature = decltype(GlobalDataFor<Feature, FeaturesGlobalData>());

  template <typename Feature>
  static inline constexpr bool contains_feature =
      tuple_contians_v<FeaturesType, std::vector<Feature>> or
      not std::is_same_v<nullptr_t, GlobalDataForFeature<Feature>>;

  template <typename DAG>
  class ViewBase
      : public std::conditional_t<
            DAG::is_mutable, typename FeatureTraits<Features, ViewType<DAG>>::Writer,
            typename FeatureTraits<Features, ViewType<DAG>>::Reader>... {};

  template <typename Feature>
  const Feature& GetFeatureAt(Id id) const;

  template <typename Feature>
  void SetFeatureAt(Id id, Feature&& feature);

  template <typename Feature>
  auto& GetFeatureGlobalData();

  template <typename Feature>
  const auto& GetFeatureGlobalData() const;

  template <size_t I = 0>
  void EnsureFeaturesSize(size_t size);

 private:
  DAG_VIEW_FRIENDS;
  [[no_unique_address]] std::tuple<std::vector<Features>...> features_;
  [[no_unique_address]] FeaturesGlobalData features_global_data_;
};