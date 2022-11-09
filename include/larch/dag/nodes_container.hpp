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

  template <typename Feature>
  Feature& GetFeatureAt(NodeId id);

  template <typename Feature>
  const Feature& GetFeatureAt(NodeId id) const;

  Storage& AddNode(NodeId id);
  void InitializeNodes(size_t nodes_count);
  auto View();
  auto View() const;
  Storage& NodeAt(NodeId id);
  const Storage& NodeAt(NodeId id) const;
  size_t Count() const;

 private:
  DAG_VIEW_FRIENDS;
  std::vector<Storage> nodes_;
  [[no_unique_address]] std::tuple<std::vector<Features>...> features_;
};