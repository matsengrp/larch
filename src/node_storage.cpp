#include "dag.hpp"

const std::vector<EdgeId>& NodeStorage::GetParents() const { return parents_; }

const std::vector<std::vector<EdgeId>>& NodeStorage::GetClades() const {
  return clades_;
}

void NodeStorage::ClearConnections() {
  parents_.clear();
  clades_.clear();
}

void NodeStorage::AddEdge(CladeIdx clade, EdgeId id, bool this_node_is_parent) {
  if (this_node_is_parent) {
    GetOrInsert(clades_, clade).push_back(id);
  } else {
    parents_.push_back(id);
  }
}

const std::optional<std::string>& NodeStorage::GetSampleId() const {
  return sample_id_;
}

void NodeStorage::SetSampleId(std::optional<std::string>&& sample_id) {
  sample_id_ = std::forward<std::optional<std::string>>(sample_id);
}

void NodeStorage::RemoveEdge(Edge edge, bool this_node_is_parent) {
  if (this_node_is_parent) {
    auto& clade = clades_.at(edge.GetClade().value);
    auto i = std::find(clade.begin(), clade.end(), edge);
    if (i != clade.end()) {
      clade.erase(i);
    }
  } else {
    auto i = std::find(parents_.begin(), parents_.end(), edge);
    if (i != parents_.end()) {
      parents_.erase(i);
    }
  }
}
