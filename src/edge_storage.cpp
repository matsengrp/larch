#include "edge_storage.hpp"

NodeId EdgeStorage::GetParent() const { return parent_; }

NodeId EdgeStorage::GetChild() const { return child_; }

CladeIdx EdgeStorage::GetClade() const { return clade_; }

void EdgeStorage::Set(NodeId parent, NodeId child, CladeIdx clade) {
  parent_ = parent;
  child_ = child;
  clade_ = clade;
}

void EdgeStorage::Set(NodeId parent, NodeId child) {
  parent_ = parent;
  child_ = child;
}
