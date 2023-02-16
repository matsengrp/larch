template <typename CRTP, typename Tag>
bool FeatureConstView<MATConversion, CRTP, Tag>::HaveMATNode() const {
  size_t id = GetMATNodeId();
  auto& node = static_cast<const CRTP&>(*this);
  return node.GetDAG().GetMAT().get_node(id) != nullptr;
}

template <typename CRTP, typename Tag>
const MAT::Node& FeatureConstView<MATConversion, CRTP, Tag>::GetMATNode() const {
  size_t id = GetMATNodeId();
  auto& node = static_cast<const CRTP&>(*this);
  auto* result = node.GetDAG().GetMAT().get_node(id);
  Assert(result != nullptr);
  return *result;
}

template <typename CRTP, typename Tag>
size_t FeatureConstView<MATConversion, CRTP, Tag>::GetMATNodeId() const {
  size_t id = GetFeatureStorage(this).mat_node_id_;
  Assert(id != NoId);
  return id;
}

template <typename CRTP, typename Tag>
MAT::Node& FeatureMutableView<MATConversion, CRTP, Tag>::GetMutableMATNode() const {
  auto& node = static_cast<const CRTP&>(*this);
  size_t id = node.GetMATNodeId();
  auto* result = node.GetDAG().GetMutableMAT().get_node(id);
  Assert(result != nullptr);
  return *result;
}

template <typename CRTP, typename Tag>
void FeatureMutableView<MATConversion, CRTP, Tag>::SetMATNodeId(size_t id) const {
  GetFeatureStorage(this).mat_node_id_ = id;
}

template <typename CRTP>
const MAT::Tree& ExtraFeatureConstView<MATConversion, CRTP>::GetMAT() const {
  auto& dag = static_cast<const CRTP&>(*this);
  return dag.template GetFeatureExtraStorage<NodeId, MATConversion>().mat_tree_;
}

template <typename CRTP>
MAT::Tree& ExtraFeatureMutableView<MATConversion, CRTP>::GetMutableMAT() const {
  auto& dag = static_cast<const CRTP&>(*this);
  return dag.template GetFeatureExtraStorage<NodeId, MATConversion>().mat_tree_;
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
      return 8;  // NOLINT
    default:
      Fail("Invalid base");
  };
}

inline auto mutations_view(MAT::Node* node) {
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

}  // namespace

template <typename CRTP>
MAT::Tree& ExtraFeatureMutableView<MATConversion, CRTP>::BuildMAT() const {
  auto& dag = static_cast<const CRTP&>(*this);
  dag.AssertUA();
  auto& tree = dag.template GetFeatureExtraStorage<NodeId, MATConversion>().mat_tree_;

  auto root_node = dag.GetRoot().GetFirstChild().GetChild();
  root_node.SetMATNodeId(root_node.GetId().value);
  // NOLINTNEXTLINE(cppcoreguidelines-owning-memory)
  auto* mat_root_node = new MAT::Node(root_node.GetId().value);

  const auto& tree_root_mutations = dag.GetRoot().GetFirstChild().GetEdgeMutations();
  mat_root_node->mutations.reserve(tree_root_mutations.size());
  for (auto [pos, muts] : tree_root_mutations) {
    Assert(pos.value != NoId);
    MAT::Mutation mat_mut("ref", static_cast<int>(pos.value),
                          EncodeBaseMAT(muts.second), EncodeBaseMAT(muts.first),
                          EncodeBaseMAT(muts.second));
    mat_root_node->mutations.push_back(mat_mut);
  }

  tree.root = mat_root_node;
  tree.register_node_serial(mat_root_node);
  BuildHelper(root_node, mat_root_node, tree);
  return tree;
}

template <typename CRTP>
void ExtraFeatureMutableView<MATConversion, CRTP>::BuildFromMAT(
    MAT::Tree&& mat, std::string_view reference_sequence) const {
  auto& dag = static_cast<const CRTP&>(*this);
  Assert(dag.IsEmpty());
  auto& tree = dag.template GetFeatureExtraStorage<NodeId, MATConversion>().mat_tree_;
  tree = std::move(mat);
  dag.SetReferenceSequence(reference_sequence);
  auto root_node = dag.AppendNode();
  root_node.SetMATNodeId(tree.root->node_id);
  BuildHelper(tree.root, root_node, dag);
  dag.BuildConnections();
  dag.AddUA(EdgeMutations{mutations_view(tree.root)});
}

template <typename CRTP>
template <typename Node>
void ExtraFeatureMutableView<MATConversion, CRTP>::BuildHelper(Node dag_node,
                                                               MAT::Node* mat_par_node,
                                                               MAT::Tree& new_tree) {
  mat_par_node->children.reserve(dag_node.GetCladesCount());
  for (auto clade : dag_node.GetClades()) {
    Assert(clade.size() == 1);
    auto edge = *clade.begin();
    const auto& mutations = edge.GetEdgeMutations();
    // NOLINTNEXTLINE(cppcoreguidelines-owning-memory)
    size_t node_id = edge.GetChild().GetId().value;
    edge.GetChild().SetMATNodeId(node_id);
    auto* node = new MAT::Node(node_id);
    new_tree.register_node_serial(node);
    node->mutations.reserve(mutations.size());
    for (auto [pos, muts] : mutations) {
      Assert(pos.value != NoId);
      MAT::Mutation mat_mut("ref", static_cast<int>(pos.value),
                            EncodeBaseMAT(muts.second), EncodeBaseMAT(muts.first),
                            EncodeBaseMAT(muts.second));
      node->mutations.push_back(mat_mut);
    }
    node->parent = mat_par_node;
    mat_par_node->children.push_back(node);
    BuildHelper(edge.GetChild(), node, new_tree);
  }
}

template <typename CRTP>
template <typename Node, typename DAG>
void ExtraFeatureMutableView<MATConversion, CRTP>::BuildHelper(MAT::Node* par_node,
                                                               Node node, DAG dag) {
  for (size_t clade_idx = 0; clade_idx < par_node->children.size(); clade_idx++) {
    MAT::Node* mat_child = par_node->children[clade_idx];
    auto child_node = dag.AppendNode();
    child_node.SetMATNodeId(mat_child->node_id);
    auto child_edge = dag.AppendEdge(node, child_node, CladeIdx{clade_idx});
    child_edge.SetEdgeMutations({mutations_view(mat_child)});
    BuildHelper(mat_child, child_node, dag);
  }
}