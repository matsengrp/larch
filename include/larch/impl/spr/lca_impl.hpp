template <typename Node>
LCA FindLCA(Node n0, Node n1) {
  auto dag = n0.GetDAG();
  Assert(dag.IsTree());
  Assert(n0.GetId() != n1.GetId());

  LCA result;
  std::set<NodeId> path0, path1;
  NodeId current_node_id0 = n0, current_node_id1 = n1;
  path0.insert(current_node_id0);
  path1.insert(current_node_id1);

  while (true) {
    if (path1.find(current_node_id0) != path1.end()) {
      result.lca = current_node_id0;
      while (dag.Get(result.path1.back()).GetParentId() != current_node_id0) {
        result.path1.pop_back();
      }
      break;
    }

    if (path0.find(current_node_id1) != path0.end()) {
      result.lca = current_node_id1;
      while (dag.Get(result.path0.back()).GetParentId() != current_node_id1) {
        result.path0.pop_back();
      }
      break;
    }

    Node current_node0 = dag.Get(current_node_id0);
    Node current_node1 = dag.Get(current_node_id1);
    if (current_node0.IsRoot() and current_node1.IsRoot()) {
      break;
    }

    if (not current_node0.IsRoot()) {
      auto parent_edge = current_node0.GetSingleParent();
      result.path0.push_back(parent_edge);
      current_node_id0 = parent_edge.GetParent();
      path0.insert(current_node_id0);
    }

    if (not current_node1.IsRoot()) {
      auto parent_edge = current_node1.GetSingleParent();
      result.path1.push_back(parent_edge);
      current_node_id1 = parent_edge.GetParent();
      path1.insert(current_node_id1);
    }
  }

  return result;
}