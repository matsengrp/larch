// Start at HypotheticalTree::GetFragment() to see how this all fits together


// updates the map in a compact genome with the supplied map. If a site from
// the new map isn't already present, adds that site. If it is already present,
// the target base is updated to the one from the new map.
CompactGenome CompactGenome::ApplyChanges(std::map<size_t, char>);

// Gets the base at provided site in the compact genome. If no mutation is
// recorded relative to the reference at that site, the reference base for that
// site is returned.
void CompactGenome::GetBase(size_t site);

// This is totally fake, just an imaginary class to describe the interface that
// I use below to describe how this could work.
class HypotheticalTreeNode {
  // GetChildren and GetParent get the nodes that would be neighbors if the SPR
  // move were applied.
  std::vector<HypotheticalTreeNode> GetChildren();
  HypotheticalTreeNode GetParent();

  CompactGenome GetOldCompactGenome() {
    if (IsNewNode()) {
      return tree_.GetTarget().GetOldCompactGenome();
    } else {
      return dag_node_->GetCompactGenome();
    }
  }

  // returns whether this node corresponds to the MAT tree root
  // (equivalent to the single child of the UA node in the sample_dag_.)
  bool IsRoot();

  // returns whether this node corresponds to the source node of the SPR move
  bool IsSourceNode();

  bool IsNewNode() {
    // dag_node_ and mat_node_ have a value if and only if
    // this node existed before the SPR move.
    return dag_node_;
  }

  // returns true if this node is an ancestor of the
  // source-target-LCA (the node returned by HypotheticalTree::GetLCA()),
  // which disqualifies a node from being a nonroot anchor
  // node (this method is used in IsNonrootAnchorNode)
  bool IsLCAAncestor();

  // This only works correctly once ComputeNewCompactGenome has been called on
  // the node, if it needs to be.
  CompactGenome GetNewCompactGenome() {
    if (new_compact_genome) {
      return new_compact_genome;
    } else {
      return GetOldCompactGenome();
    }
  }

  // return a set of site indices at which there are fitch set changes
  std::set<size_t> GetSitesWithChangedFitchSets() {
    std::set<size_t> sites_with_changes;
    std::optional<std::map<size_t, MAT::Mutation_Count_Change&>&> fitch_set_map = GetFitchSetParts().second;
    if (fitch_set_map) {
      for (auto map_pair : *fitch_set_map) {
        sites_with_changes.insert(map_pair.first);
      }
    }
    return sites_with_changes;
  }

  std::pair<MAT::Mutations_Collection&, std::optional<std::map<size_t, MAT::Mutation_Count_Change&>&>> GetFitchSetParts() {
    if (IsTargetNode()){
      //then fitch sets can't have changed, but the fitch sets recorded in
      //tree_'s changed fitch set map relative to this node are meant for this
      //node's new parent!
      return {mat_node_.mutations, {}};
    } else if (IsNewNode()) {
      // Then fitch set changes are relative to the target node
      return {tree_.GetTarget().mat_node_.mutations, tree_.changed_fitch_set_map.find(tree_.GetTarget().mat_node_ *)};
    } else {
      return {mat_node_.mutations, tree_.changed_fitch_set_map.find(mat_node*)};
    }
  }

  // get the (possibly modified) fitch set at this node at the provided site.
  FitchSet GetFitchSet(size_t site){
    auto [old_fitch_sets, changes] = GetFitchSetParts();
    if (old_fitch_sets.find(site) == old_fitch_sets.mutations.end()) {
      // if no fitch set is recorded on the corresponding MAT node, we can use
      // a singleton set containing the base at this site in the parent compact genome
      // TODO: Does this work if IsSourceNode()?
      return FitchSet({GetParent().GetNewCompactGenome().GetBase(site)});
    } else if (changes) {
      return FitchSet((old_fitch_sets.find(site).get_all_major_allele() & (~changes->get_decremented())) | changes->get_incremented());
    } else {
      return FitchSet(old_fitch_sets.find(site).get_all_major_allele());
    }
  }

  // Most of the time this can just return the parent node's
  // changed_base_sites. However, it's different if the node in question is the
  // source node!
  std::set<size_t> GetParentChangedBaseSites(){
    if (self.IsSourceNode()) {
      // this node used to be below the old source parent, so we need
      // differences between the new parent's new cg, and the old cg of the old
      // parent of the source node.
      CompactGenome old_parent_cg = tree_.GetOldSourceParent().GetOldCompactGenome();
      // Parent must be the new node:
      CompactGenome new_parent_cg = GetParent().GetNewCompactGenome();
      // Imaginary method DifferingSites returns sites at which new_parent_cg
      // and old_parent_cg don't have the same base.
      return old_parent_cg.DifferingSites(new_parent_cg);
    } else {
      return GetParent()._changed_base_sites;
    }
  }

  void ComputeNewCompactGenome() {
    std::set<size_t>& changed_parent_sites = GetSitesWithChangedFitchSets();
    CompactGenome& old_cg = GetOldCompactGenome();
    std::map<size_t, char> cg_changes;
    if IsRoot() {
      // If this node is the root node of the tree, we don't have to worry
      // about parent bases, we just
      // choose a base from each changed fitch set, preferring the base that was
      // already in that site before the SPR move
      std::set<size_t> focus_sites = GetSitesWithChangedFitchSets();
      for (auto site : focus_sites) {
        FitchSet site_fitch_set = GetFitchSet(site);
        char oldbase = old_cg.GetBase(site);
        if (not site_fitch_set.find(oldbase)) {
          cg_changes.insert({site, site_fitch_set.at(0)});
          changed_base_sites.insert(site);
        }
      }
    } else {
      // If this node is not the root node, we do need to prefer the base
      // from the new compact genome of the parent node. If it's not in the
      // fitch set, then we're free to choose any base in the fitch set, so we
      // choose the base in this node's old compact genome when that base is in
      // the fitch set.
      std::set<size_t> focus_sites = GetSitesWithChangedFitchSets();
      focus_sites.merge(GetParentChangedBaseSites());
      for (auto site : focus_sites) {
        FitchSet site_fitch_set = GetFitchSet(site);
        char oldbase = old_cg.GetBase(site);
        char parent_base = GetParent().GetCompactGenome().GetBase(site);
        if (site_fitch_set.find(parent_base)) {
          if (oldbase != parent_base) {
            cg_changes.insert({site, parent_base});
            changed_base_sites.insert(site);
          }
        } else if (not site_fitch_set.find(oldbase)) {
          cg_changes.insert({site, site_fitch_set.at(0)});
          changed_base_sites.insert(site);
        }
      }
    }
    return old_cg.ApplyChanges(cg_changes);
  }

  void PreorderComputeCompactGenome(std::vector<HypotheticalTreeNode>& fragment_vector) {
    ComputeNewCompactGenome();
    fragment_vector.push_back(this);
    // If we've reached an anchor node, there's no need to continue down this
    // branch.
    if (not self.IsNonrootAnchorNode()) {
      for (auto child : GetChildNodes()) {
        child.PreorderComputeCompactGenome();
      }
    }
  }

  // A tree fragment has two kinds of anchor nodes: One root anchor node, which
  // is always the parent of HypotheticalTree.GetOldestChangedNode(), and the
  // rest are what I'm calling here "nonroot anchor nodes", whose descendants
  // in the hypothetical tree are guaranteed to also be unchanged from before
  // the SPR move. This method identifies the second kind.
  bool IsNonrootAnchorNode() {
    if (not IsLCAAncestor()) {
      return (dag_node_.GetCompactGenome() == GetNewCompactGenome() and dag_node_.GetLeafSet() == GetNewLeafSet());
    } else {
      return false;
    }
  }

  std::optional<CompactGenome> new_compact_genome;
  std::set<size_t> changed_base_sites;
  private:
    HypotheticalTree tree_;
    std::optional<Node&> dag_node_;
    std::optional<MAT::Node&> mat_node_;
    MAT::changed_major_allele_set& fitch_set_changes;
}

// HypotheticalTree is the wrapper for an unmodified sampled tree (DAG) and the
// corresponding unmodified MAT, and for a proposed SPR move.
class HypotheticalTree {

  HypotheticalTree(const MADAG& sample_tree, const MAT& sample_mat, const Profitable_Moves& move, MAT::Node_With_Major_Allele_Set_Change& nodes_with_major_allele_set_change)
    : sample_dag_{sample_tree},
      sample_mat_{sample_mat},
      move_{move} {
        for (auto node_with_allele_set_change : nodes_with_major_allele_set_change) {
          std::map<size_t, MAT::Mutation_Count_Change> node_map;
          for (auto mutation_count_change : node_with_allele_set_change.major_allele_set_change) {
            node_map.push_back({mutation_count_change.get_position(), mutation_count_change})
          }
          changed_fitch_set_map_.insert({node_with_allele_set_change.node, node_with_allele_set_change.major_allele_set_change})
        }
      }

  // Get the LCA of source and target nodes
  HypotheticalTreeNode GetLCA();
  // These return the HypotheticalTreeNodes corresponding to source and target
  // nodes (they're siblings in the hypothetical tree)
  HypotheticalTreeNode GetSource();
  HypotheticalTreeNode GetTarget();
  
  // Returns the HypotheticalTreeNode that used to be the parent of source
  // before the SPR move. TODO: This node may (but need not be) unifurcating
  // after the SPR move, in which case the issue description says it should not
  // appear as a neighbor of any node in the hypothetical tree (it's omitted
  // entirely). However, we will still need to access its data, such as its old
  // compact genome, so this method should return it. the todo is to think
  // through all the consequences of this, and make sure it doesn't cause any
  // problems that we skip it sometimes in the tree traversal.
  HypotheticalTreeNode GetOldSourceParent();

  // returns LCA of source and target, or the earliest node with fitch set
  // changes, whichever is higher in the tree.
  HypotheticalTreeNode GetOldestChangedNode();

  std::vector<HypotheticalTreeNode> GetFragment() {
    std::vector<HypotheticalTreeNode> result;
    HypotheticalTreeNode& oldest_changed = GetOldestChangedNode();
    if (oldest_changed.IsRoot()) {
      // we need to add the UA node as the root anchor node of the fragment,
      // somehow
    } else {
      result.push_back(oldest_changed.GetParent());
    }
    oldest_changed.PreorderComputeCompactGenome(result);
    return result;
  }

  private:
    const MADAG& sample_dag_;
    const MAT& sample_mat_;
    const Profitable_Moves& move;
    std::map<MAT::Node*, std::map<size_t, MAT::Mutation_Count_Change&>> changed_fitch_set_map_;
}

// // From
// // usher/src/matOptimize/Profitable_MovesEnumerators/Profitable_Moves_Enumerators.hpp
// // for reference here are some relevant data structures...
//
// std::vector<Node_With_Major_Allele_Set_Change>& node_with_major_allele_set_change;
// Node_With_Major_Allele_Set_Change {
//   MAT::Node* node;
//   Mutation_Count_Change_Collection major_allele_set_change;
// }
// 
// typedef std::vector<Mutation_Count_Change> Mutation_Count_Change_Collection;
// 
//
// // Also, usher/src/matOptimize/mutation_annotated_tree.hpp contains
// definition for MAT::Node, whose mutations attribute holds
// a MAT::Mutations_Collection object, which is a wrapper for a vector of
// MAT::Mutation objects.
