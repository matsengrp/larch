#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

template <typename Storage, typename... Features>
class DefaultEdgesContainer
    : public DefaultContainerBase<EdgeId, EdgeView, Storage, Features...> {
 public:
  Storage& AddEdge(EdgeId id);
  auto View();
  auto View() const;
  size_t Count() const;
  Storage& EdgeAt(EdgeId id);
  const Storage& EdgeAt(EdgeId id) const;

 private:
  DAG_VIEW_FRIENDS;
  std::vector<Storage> edges_;
};