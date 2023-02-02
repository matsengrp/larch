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
  bool IsRoot()

  bool IsNewNode() {
    // dag_node_ and mat_node_ have a value if and only if
    // this node existed before the SPR move.
    return dag_node_;
  }

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
  std::set<size_t> GetSitesWithChangedFitchSets();

  // get the (possibly modified) fitch set at this node at the provided site.
  FitchSet GetFitchSet(size_t site){
    // if no fitch set is recorded on the corresponding MAT node, we can use
    // a singleton set containing the base at this site in the parent compact genome
    if (mat_node_.fitch_sets.at(site).empty()) {
      return FitchSet({GetParent().GetNewCompactGenome().GetBase(site)})
    } else {
      return FitchSet((mat_node_.fitch_sets.at(site) & (~fitch_set_changes.decr_alleles)) | fitch_set_changes.incr_alleles)
    }
  }

  void ComputeNewCompactGenome() {
    std::vector<size_t>& changed_parent_sites;
    CompactGenome& old_cg = GetOldCompactGenome();
    std::map<size_t, char> cg_changes;
    if IsRoot() {
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
      std::set<size_t> focus_sites = GetSitesWithChangedFitchSets();
      focus_sites.merge(GetParent().changed_base_sites);
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

  bool IsAnchorNode();

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

  HypotheticalTree(const MADAG& sample_tree, const MAT& sample_mat, const Profitable_Moves& move)
    : sample_dag_{sample_tree},
      sample_mat_{sample_mat},
      move_{move} {}

  // Get the LCA of source and target nodes
  HypotheticalTreeNode GetLCA();
  // These return the HypotheticalTreeNodes corresponding to source and target
  // nodes (they're siblings in the hypothetical tree)
  HypotheticalTreeNode GetSource();
  HypotheticalTreeNode GetTarget();
  
  // Returns the HypotheticalTreeNode that used to be the parent of source
  // before the SPR move.
  HypotheticalTreeNode GetOldSourceParent();

  // returns LCA of source and target, or the earliest node with fitch set
  // changes, whichever is higher in the tree.
  HypotheticalTreeNode GetOldestChangedNode();

  private:
    const MADAG& sample_dag_;
    const MAT& sample_mat_;
    const Profitable_Moves& move;
}


