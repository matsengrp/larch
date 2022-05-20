#include <iostream>  //XXX

Merge::Merge(const std::string_view& refseq,
             std::vector<std::reference_wrapper<const HistoryDAG>>&& trees,
             const std::vector<TreeLabels>& labels)
    : refseq_{refseq}, trees_{std::move(trees)}, labels_{labels} {
  assert(labels_.size() >= trees_.size());
  for (size_t tree_idx = 0; tree_idx < trees_.size(); ++tree_idx) {
    assert(labels_.at(tree_idx).compact_genomes.size() >=
           trees_.at(tree_idx).get().GetNodes().size());
    assert(labels_.at(tree_idx).leaf_sets.size() >=
           trees_.at(tree_idx).get().GetNodes().size());
  }
}

void Merge::Run() {
  std::cout << "Merging trees ";
  for (size_t tree_idx = 0; tree_idx < trees_.size(); ++tree_idx) {
    std::cout << "." << std::flush;
    for (Edge edge : trees_.at(tree_idx).get().GetEdges()) {
      MakeResultEdge(tree_idx, edge.GetId());
    }
  }
  std::cout << " done."
            << "\n";
  result_.BuildConnections();
}

NodeLabel Merge::GetNodeLabel(size_t tree_idx, NodeId node_id) {
  return labels_.at(tree_idx).GetLabel(node_id);
}

NodeId Merge::GetResultNode(NodeLabel label) {
  auto i = result_nodes_.find(label);
  NodeId id;
  if (i == result_nodes_.end()) {
    id = result_.AddNode({result_.GetNodes().size()}).GetId();
    result_nodes_.emplace_hint(i, label, id);
  } else {
    id = i->second;
  }
  return id;
}

void Merge::MakeResultEdge(size_t tree_idx, EdgeId edge_id) {
  Edge edge = trees_.at(tree_idx).get().GetEdge(edge_id);
  NodeId parent = edge.GetParent().GetId(), child = edge.GetChild().GetId();
  EdgeLabel label = {
      GetNodeLabel(tree_idx, parent), GetNodeLabel(tree_idx, child), {0}};
  auto i = result_edges_.find(label);
  if (i == result_edges_.end()) {
    NodeId result_parent = GetResultNode(std::get<0>(label));
    NodeId result_child = GetResultNode(std::get<1>(label));
    result_.AddEdge({result_.GetEdges().size()}, result_parent, result_child,
                    std::get<2>(label));
    result_edges_.emplace_hint(i, std::move(label));
  }
}

HistoryDAG& Merge::GetResult() { return result_; }

const HistoryDAG& Merge::GetResult() const { return result_; }

TreeLabels GetLabels(const HistoryDAG& tree, const std::string_view& refseq,
                     const std::vector<CompactGenome>& mutations) {
  TreeLabels result;
  result.compact_genomes.resize(tree.GetNodes().size());
  for (auto iter : tree.TraversePreOrder()) {
    if (iter.IsRoot()) {
      continue;
    }
    const CompactGenome& muts = mutations.at(iter.GetEdge().GetId().value);
    const CompactGenome& parent_cgs =
        result.compact_genomes.at(iter.GetEdge().GetParent().GetId().value);
    CompactGenome& cgs = result.compact_genomes.at(iter.GetNode().GetId().value);
    cgs = parent_cgs;
    for (auto [pos, base] : muts) {
      if (base != refseq.at(pos - 1)) {
        cgs[pos] = base;
      } else {
        cgs.erase(pos);
      }
    }
  }

  result.leaf_sets.resize(tree.GetNodes().size());
  for (Node node : tree.TraversePostOrder()) {
    if (node.IsLeaf()) {
      continue;
    }
    LeafSet& leaf_set = result.leaf_sets.at(node.GetId().value);
    for (auto clade : node.GetClades()) {
      std::unordered_set<CompactGenomePointer> clade_leafs;
      for (Node child : clade | ranges::views::transform(Transform::GetChild)) {
        const LeafSet& child_leaf_set = result.leaf_sets.at(child.GetId().value);
        if (child.IsLeaf()) {
          clade_leafs.insert({&result.compact_genomes.at(child.GetId().value)});
        } else {
          for (auto& child_leafs : result.leaf_sets.at(child.GetId().value)) {
            clade_leafs.insert(child_leafs.begin(), child_leafs.end());
          }
        }
      }
      leaf_set.push_back(clade_leafs);
    }
  }

  return result;
}
