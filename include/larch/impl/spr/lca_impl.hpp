template <typename Node>
LCA FindLCA(Node n0, Node n1) {
  static_assert(Node::role == Role::View);
  static_assert(Node::component == Component::Node);
  Assert(std::addressof(n0.GetDAG().GetStorage()) ==
         std::addressof(n1.GetDAG().GetStorage()));
  auto dag = n0.GetDAG();
  Assert(dag.IsTree());
  Assert(n0.GetId() != n1.GetId());

  ContiguousSet<NodeId> visited;
  LCA result;

  auto visit = [dag, &visited, &result](EdgeId edge,
                                        std::vector<EdgeId>& path) -> EdgeId {
    path.push_back(edge);
    auto parent = dag.Get(edge).GetParent();
    if (visited.insert(parent).second) {
      if (not parent.IsUA()) {
        return parent.GetSingleParent();
      } else {
        return {NoId};
      }
    } else {
      result.lca = parent;
      return {NoId};
    }
  };

  EdgeId e0 = n0.IsUA() ? EdgeId{NoId} : n0.GetSingleParent();
  EdgeId e1 = n1.IsUA() ? EdgeId{NoId} : n1.GetSingleParent();

  visited.insert(n0);
  visited.insert(n1);

  while (e0.value != NoId or e1.value != NoId) {
    if (e0.value != NoId) {
      e0 = visit(e0, result.path0);
    }
    if (e1.value != NoId) {
      e1 = visit(e1, result.path1);
    }
  }

  Assert(result.lca.value != NoId);

  return result;
}
