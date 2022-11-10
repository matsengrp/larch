#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

template <typename Storage, typename... Features>
class DefaultNodesContainer
    : public DefaultContainerBase<NodeId, NodeView, Storage, Features...> {
 public:
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
};