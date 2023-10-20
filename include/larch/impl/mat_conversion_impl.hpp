template <typename CRTP, typename Tag>
bool FeatureConstView<MATConversion, CRTP, Tag>::HaveMATNode() const {
  MATNodePtr ptr = GetFeatureStorage(this).mat_node_ptr_;
  return ptr != nullptr;
}

template <typename CRTP, typename Tag>
bool FeatureConstView<MATConversion, CRTP, Tag>::IsCondensedInMAT() const {
  return GetFeatureStorage(this).is_condensed_in_mat_;
}

template <typename CRTP, typename Tag>
MATNodePtr FeatureConstView<MATConversion, CRTP, Tag>::GetMATNode() const {
  MATNodePtr ptr = GetFeatureStorage(this).mat_node_ptr_;
  Assert(ptr != nullptr);
  return ptr;
}

template <typename CRTP, typename Tag>
MATNodePtr FeatureMutableView<MATConversion, CRTP, Tag>::GetMutableMATNode() const {
  MATNodePtr ptr = GetFeatureStorage(this).mat_node_ptr_;
  Assert(ptr != nullptr);
  return ptr;
}

template <typename CRTP, typename Tag>
void FeatureMutableView<MATConversion, CRTP, Tag>::SetMATNode(MATNodePtr ptr) const {
  GetFeatureStorage(this).mat_node_ptr_ = ptr;
  auto& node = static_cast<const CRTP&>(*this);
  node.GetDAG()
      .template GetFeatureExtraStorage<Component::Node, MATConversion>()
      .reverse_map_.insert({ptr, node.GetId()});
}

template <typename CRTP, typename Tag>
void FeatureMutableView<MATConversion, CRTP, Tag>::SetUncondensedMATNode(
    MATNodePtr ptr) const {
  GetFeatureStorage(this).mat_node_ptr_ = ptr;
  auto& node = static_cast<const CRTP&>(*this);

  if (GetFeatureStorage(this).is_condensed_in_mat_) {
    node.GetDAG()
        .template GetFeatureExtraStorage<Component::Node, MATConversion>()
        .uncondensed_reverse_map_.at(ptr)
        .push_back(node.GetId());
  } else {
    GetFeatureStorage(this).is_condensed_in_mat_ = true;
    std::vector<NodeId> nodevec;
    nodevec.push_back(node.GetId());
    node.GetDAG()
        .template GetFeatureExtraStorage<Component::Node, MATConversion>()
        .uncondensed_reverse_map_.insert({ptr, nodevec});
  }
}

template <typename CRTP>
const MAT::Tree& ExtraFeatureConstView<MATConversion, CRTP>::GetMAT() const {
  auto& dag = static_cast<const CRTP&>(*this);
  return *dag.template GetFeatureExtraStorage<Component::Node, MATConversion>()
              .mat_tree_;
}

template <typename CRTP>
auto ExtraFeatureConstView<MATConversion, CRTP>::GetNodeFromMAT(MATNodePtr node) const {
  auto& dag = static_cast<const CRTP&>(*this);
  NodeId id = dag.template GetFeatureExtraStorage<Component::Node, MATConversion>()
                  .reverse_map_.at(node);
  return dag.Get(id);
}

template <typename CRTP>
auto ExtraFeatureConstView<MATConversion, CRTP>::GetUncondensedNodeFromMAT(
    MATNodePtr node) const {
  auto& dag = static_cast<const CRTP&>(*this);
  if (dag.Get(dag.template GetFeatureExtraStorage<Component::Node, MATConversion>()
                  .reverse_map_.at(node))
          .IsCondensedInMAT()) {
    return dag.template GetFeatureExtraStorage<Component::Node, MATConversion>()
        .uncondensed_reverse_map_.at(node);
  }
  std::vector<NodeId> to_ret;
  to_ret.push_back(dag.template GetFeatureExtraStorage<Component::Node, MATConversion>()
                       .reverse_map_.at(node));
  return to_ret;
}

template <typename CRTP>
MAT::Tree& ExtraFeatureMutableView<MATConversion, CRTP>::GetMutableMAT() const {
  auto& dag = static_cast<const CRTP&>(*this);
  auto* result =
      dag.template GetFeatureExtraStorage<Component::Node, MATConversion>().mat_tree_;
  Assert(result != nullptr);
  return *result;
}

template <typename CRTP>
auto ExtraFeatureMutableView<MATConversion, CRTP>::GetMutableNodeFromMAT(
    MATNodePtr node) const {
  auto& dag = static_cast<const CRTP&>(*this);
  NodeId id = dag.template GetFeatureExtraStorage<Component::Node, MATConversion>()
                  .reverse_map_.at(node);
  return dag.Get(id);
}

namespace {

static inline uint8_t EncodeBaseMAT(char base) {
  switch (base) {
    case 'A':
      return 1;
    case 'C':
      return 2;
    case 'G':
      return 4;
    case 'T':
      return 8;
    case 'N':
      return 1 + 2 + 4 + 8;  // NOLINT
    default:
      Fail("Invalid base");
  };
}

inline auto mutations_view(MATNodePtr node) {
  return node->mutations |
         ranges::views::transform(
             [](const MAT::Mutation& mut)
                 -> std::pair<MutationPosition, std::pair<char, char>> {
               static const std::array<char, 4> decode = {'A', 'C', 'G', 'T'};
               return {{static_cast<size_t>(mut.get_position())},
                       {decode.at(one_hot_to_two_bit(mut.get_par_one_hot())),
                        decode.at(one_hot_to_two_bit(mut.get_mut_one_hot()))}};
             });
}

static inline void fill_static_reference_sequence(std::string_view dag_ref) {
  static std::mutex static_ref_seq_mutex;
  std::lock_guard lock{static_ref_seq_mutex};
  MAT::Mutation::refs.resize(dag_ref.size() + 1);
  for (size_t ref_idx = 0; ref_idx < dag_ref.size(); ref_idx++) {
    MAT::Mutation::refs[ref_idx + 1] = EncodeBaseMAT(dag_ref[ref_idx]);
  }
}

}  // namespace

template <typename CRTP>
void ExtraFeatureMutableView<MATConversion, CRTP>::BuildMAT(MAT::Tree& tree) const {
  auto& dag = static_cast<const CRTP&>(*this);
  dag.AssertUA();
  fill_static_reference_sequence(dag.GetReferenceSequence());
  dag.template GetFeatureExtraStorage<Component::Node, MATConversion>().mat_tree_ =
      std::addressof(tree);

  auto root_node = dag.GetRoot().GetFirstChild().GetChild();
  // NOLINTNEXTLINE(cppcoreguidelines-owning-memory)
  MATNodePtr mat_root_node = new MAT::Node(root_node.GetId().value);
  root_node.SetMATNode(mat_root_node);

  const auto& tree_root_mutations = dag.GetRoot().GetFirstChild().GetEdgeMutations();
  mat_root_node->mutations.reserve(tree_root_mutations.size());
  for (auto [pos, muts] : tree_root_mutations) {
    Assert(pos.value != NoId);
    MAT::Mutation mat_mut(
        "ref", static_cast<int>(pos.value), EncodeBaseMAT(muts.second.ToChar()),
        EncodeBaseMAT(muts.first.ToChar()), EncodeBaseMAT(muts.second.ToChar()));
    mat_root_node->mutations.push_back(mat_mut);
  }

  tree.root = mat_root_node;
  tree.register_node_serial(mat_root_node);
  BuildHelper(root_node, mat_root_node, tree);
}

template <typename CRTP>
void ExtraFeatureMutableView<MATConversion, CRTP>::BuildFromMAT(
    MAT::Tree& mat, std::string_view reference_sequence) const {
  auto& dag = static_cast<const CRTP&>(*this);
  Assert(dag.IsEmpty());
  dag.template GetFeatureExtraStorage<Component::Node, MATConversion>().mat_tree_ =
      std::addressof(mat);
  dag.SetReferenceSequence(reference_sequence);
  auto root_node = dag.AppendNode();
  root_node.SetMATNode(mat.root);
  BuildHelper(mat.root, root_node, dag);
  dag.BuildConnections();

  if (mat.condensed_nodes.size() > 0) {
    for (auto cn = mat.condensed_nodes.begin(); cn != mat.condensed_nodes.end(); cn++) {
      auto node_to_uncondense = mat.get_node(cn->first);
      auto parent_node = (node_to_uncondense->parent != NULL)
                             ? node_to_uncondense->parent
                             : node_to_uncondense;
      auto dag_parent_node = dag.GetNodeFromMAT(
          parent_node);  // there's only ever one parent because MATs are tree-shaped
      size_t num_samples = cn->second.size();
      if (num_samples > 0) {
        // reset the DAG node that points to the condensed node so that it
        // now points to the first node in the vector of condensed nodes.
        for (auto node : dag_parent_node.GetChildren() | Transform::GetChild()) {
          if (node.GetMATNode() == node_to_uncondense) {
            node.SetMATNode(mat.get_node(cn->second[0]));
            node.SetUncondensedMATNode(mat.get_node(cn->second[0]));
            break;
          }
        }
        // add all of the remaining condensed nodes as siblings
        CladeIdx clade_idx = {dag_parent_node.GetCladesCount()};
        for (size_t s = 1; s < num_samples; s++) {
          auto dag_child_node = dag.AppendNode();
          auto mat_child_node = mat.get_node(cn->second[s]);
          dag_child_node.SetMATNode(mat.get_node(cn->second[0]));
          dag_child_node.SetUncondensedMATNode(mat_child_node);

          auto child_edge = dag.AppendEdge(dag_parent_node, dag_child_node, clade_idx);
          child_edge.SetEdgeMutations({mutations_view(mat_child_node)});
          ++clade_idx.value;
        }
      }
    }
    dag.BuildConnections();
  }
  for (auto leaf : dag.GetLeafs()) {
    if constexpr (decltype(leaf)::template contains_feature<Deduplicate<SampleId>>) {
      auto id_iter = dag.template AsFeature<Deduplicate<SampleId>>().AddDeduplicated(
          SampleId{std::to_string(leaf.GetMATNode()->node_id)});
      leaf = id_iter.first;
    } else {
      leaf.SetSampleId(std::to_string(leaf.GetMATNode()->node_id));
    }
  }
  dag.AddUA(EdgeMutations{mutations_view(mat.root)});
}

template <typename CRTP>
template <typename Node>
void ExtraFeatureMutableView<MATConversion, CRTP>::BuildHelper(Node dag_node,
                                                               MATNodePtr mat_par_node,
                                                               MAT::Tree& new_tree) {
  mat_par_node->children.reserve(dag_node.GetCladesCount());
  for (auto clade : dag_node.GetClades()) {
    Assert(clade.size() == 1);
    auto edge = *clade.begin();
    const auto& mutations = edge.GetEdgeMutations();
    // NOLINTNEXTLINE(cppcoreguidelines-owning-memory)
    size_t node_id = edge.GetChild().GetId().value;
    auto node_name = edge.GetChild().GetSampleId().value_or(std::to_string(node_id));
    auto* node = new MAT::Node(node_id);
    edge.GetChild().SetMATNode(node);
    node->clade_annotations.resize(new_tree.get_num_annotations(), "");
    new_tree.register_node_serial(node, node_name);
    node->mutations.reserve(mutations.size());
    for (auto [pos, muts] : mutations) {
      Assert(pos.value != NoId);
      MAT::Mutation mat_mut(
          "ref", static_cast<int>(pos.value), EncodeBaseMAT(muts.second.ToChar()),
          EncodeBaseMAT(muts.first.ToChar()), EncodeBaseMAT(muts.second.ToChar()));
      node->mutations.push_back(mat_mut);
    }
    node->parent = mat_par_node;
    mat_par_node->children.push_back(node);
    BuildHelper(edge.GetChild(), node, new_tree);
  }
}

namespace TODO {
static inline std::string ToEdgeMutationsString(MATNodePtr node) {  // XXX
  static const std::array<char, 4> decode = {'A', 'C', 'G', 'T'};
  std::string result = "<";
  for (const MAT::Mutation& mut : node->mutations) {
    result += decode.at(one_hot_to_two_bit(mut.get_par_one_hot()));
    result += std::to_string(mut.get_position());
    result += decode.at(one_hot_to_two_bit(mut.get_mut_one_hot()));
    result += ", ";
  }
  return result + ">";
}
}  // namespace TODO

template <typename CRTP>
template <typename Node, typename DAG>
void ExtraFeatureMutableView<MATConversion, CRTP>::BuildHelper(MATNodePtr par_node,
                                                               Node node, DAG dag) {
  for (CladeIdx clade_idx = {0}; clade_idx.value < par_node->children.size();
       ++clade_idx.value) {
    MATNodePtr mat_child = par_node->children[clade_idx.value];
    auto child_node = dag.AppendNode();
    child_node.SetMATNode(mat_child);
    auto child_edge = dag.AppendEdge(node, child_node, clade_idx);
    child_edge.SetEdgeMutations({mutations_view(mat_child)});
    BuildHelper(mat_child, child_node, dag);
  }
}
