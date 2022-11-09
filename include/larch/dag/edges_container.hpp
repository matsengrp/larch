#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

template <typename Storage, typename... Features>
class DefaultEdgesContainer {
 public:
  using StorageType = Storage;
  using FeaturesType = std::tuple<std::vector<Features>...>;

  template <typename Feature>
  static inline constexpr bool contains_feature =
      tuple_contians_v<FeaturesType, std::vector<Feature>>;

  template <typename DAG>
  class ViewBase
      : public std::conditional_t<DAG::is_mutable,
                                  FeatureWriter<Features, EdgeView<DAG>>,
                                  FeatureReader<Features, EdgeView<DAG>>>... {};

  template <typename Feature>
  Feature& GetFeatureAt(EdgeId id);

  template <typename Feature>
  const Feature& GetFeatureAt(EdgeId id) const;

  Storage& AddEdge(EdgeId id);
  auto View();
  auto View() const;
  size_t Count() const;
  Storage& EdgeAt(EdgeId id);
  const Storage& EdgeAt(EdgeId id) const;

 private:
  DAG_VIEW_FRIENDS;
  std::vector<Storage> edges_;
  [[no_unique_address]] std::tuple<std::vector<Features>...> features_;
};