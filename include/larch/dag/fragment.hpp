#ifndef DAG_DECLARATIONS
#error "Don't include this header, use larch/dag/dag.hpp instead"
#endif

template <typename DAG>
class FragmentView;

template <typename DAG>
class FragmentStorage {
 public:
  MOVE_ONLY(FragmentStorage);
  FragmentStorage(DAG dag, std::vector<NodeId>&& nodes, std::vector<EdgeId>&& edges)
      : dag_{dag}, nodes_{std::move(nodes)}, edges_{std::move(edges)} {}
  auto View() { return FragmentView{*this}; }
  auto View() const { return FragmentView{*this}; }

 private:
  friend class FragmentView<DAG>;
  DAG dag_;
  const std::vector<NodeId> nodes_;
  const std::vector<EdgeId> edges_;
};

template <typename DAG>
class FragmentView {
 public:
  template <Component C, typename Feature>
  static const bool contains_element_feature =
      DAG::template contains_element_feature<C, Feature>;

  explicit FragmentView(FragmentStorage<DAG>& storage);

  void AssertUA() const;
  size_t GetNodesCount() const;
  size_t GetEdgesCount() const;
  auto Get(NodeId id) const;
  auto Get(EdgeId id) const;
  auto GetNodes() const;
  auto GetEdges() const;
  auto GetRoot() const;

 private:
  FragmentStorage<DAG>& storage_;
};
