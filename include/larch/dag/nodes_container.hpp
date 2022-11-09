#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

template <typename Storage, typename... Features>
class DefaultNodesContainer {
 public:
  using StorageType = Storage;
  using FeaturesType = std::tuple<std::vector<Features>...>;

  template <typename Feature>
  static inline constexpr bool contains_feature =
      tuple_contians_v<FeaturesType, std::vector<Feature>>;

  template <typename DAG>
  class ViewBase
      : public std::conditional_t<DAG::is_mutable,
                                  FeatureWriter<Features, NodeView<DAG>>,
                                  FeatureReader<Features, NodeView<DAG>>>... {};

  template <typename Feature, typename Id>
  Feature& GetFeatureAt(Id id);

  template <typename Feature, typename Id>
  const Feature& GetFeatureAt(Id id) const;

 private:
  DAG_VIEW_FRIENDS;
  std::vector<Storage> nodes_;
  [[no_unique_address]] std::tuple<std::vector<Features>...> features_;
};