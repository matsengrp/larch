#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

template <typename DAG>
class EdgeView
    : public DAG::StorageType::EdgesContainerType::StorageType::template ViewBase<DAG>,
      public DAG::StorageType::EdgesContainerType::template ViewBase<DAG> {
 public:
  constexpr static const bool is_mutable = DAG::is_mutable;
  using Node = typename DAG::Node;
  EdgeView(DAG dag, EdgeId id);
  operator EdgeView<typename DAG::Immutable>();
  operator EdgeId();
  operator CladeIdx();
  auto& GetDAG();
  EdgeId GetId();
  auto GetParent();
  auto GetChild();
  CladeIdx GetClade();
  NodeId GetParentId();
  NodeId GetChildId();
  std::pair<NodeId, NodeId> GetNodeIds();
  bool IsRoot();
  bool IsLeaf();

 private:
  DAG_FEATURE_FRIENDS;
  auto& GetStorage() const;
  template <typename Feature>
  auto& GetFeatureStorage() const;
  DAG dag_;
  EdgeId id_;
};

namespace std {
template <typename DAG>
struct tuple_size<::EdgeView<DAG>> : integral_constant<size_t, 2> {};

template <size_t Index, typename DAG>
struct tuple_element<Index, ::EdgeView<DAG>>
    : tuple_element<Index, tuple<typename DAG::Node, typename DAG::Node>> {};
}  // namespace std

template <std::size_t Index, typename DAG>
std::tuple_element_t<Index, EdgeView<DAG>> get(EdgeView<DAG> edge);