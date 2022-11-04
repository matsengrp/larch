#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

template <typename DAGType, typename... Features>
class EdgeView
    : public std::conditional_t<
          DAGType::is_mutable, FeatureWriter<Features, EdgeView<DAGType, Features...>>,
          FeatureReader<Features, EdgeView<DAGType, Features...>>>... {
 public:
  constexpr static const bool is_mutable = DAGType::is_mutable;
  using Node = typename DAGType::Node;
  EdgeView(DAGType dag, EdgeId id);
  operator EdgeView<typename DAGType::Immutable, Features...>() const;
  operator EdgeId() const;
  operator CladeIdx() const;
  auto& GetDAG() const;
  EdgeId GetId() const;
  auto GetParent() const;
  auto GetChild() const;
  CladeIdx GetClade() const;
  NodeId GetParentId() const;
  NodeId GetChildId() const;
  std::pair<NodeId, NodeId> GetNodeIds() const;
  bool IsRoot() const;
  bool IsLeaf() const;

 private:
  DAG_FEATURE_FRIENDS;
  auto& GetStorage() const;
  DAGType dag_;
  EdgeId id_;
};

namespace std {
template <typename DAGType, typename... Features>
struct tuple_size<::EdgeView<DAGType, Features...>> : integral_constant<size_t, 2> {};

template <size_t Index, typename DAGType, typename... Features>
struct tuple_element<Index, ::EdgeView<DAGType, Features...>>
    : tuple_element<Index, tuple<typename DAGType::Node, typename DAGType::Node>> {};
}  // namespace std

template <std::size_t Index, typename DAGType, typename... Features>
std::tuple_element_t<Index, EdgeView<DAGType, Features...>> get(
    EdgeView<DAGType, Features...> edge);