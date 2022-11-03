#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

template <typename Storage>
class DefaultNodesContainer {
 public:
  using StorageType = Storage;

 private:
  DAG_VIEW_FRIENDS;
  std::vector<Storage> nodes_;
};