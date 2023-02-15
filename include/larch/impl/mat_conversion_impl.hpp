template <typename CRTP, typename Tag>
const MAT::Node& FeatureConstView<MATConversion, CRTP, Tag>::GetMATNode() const {
  auto& node = static_cast<const CRTP&>(*this);
  auto& mat = node.GetDAG().GetMAT();
  size_t id = GetFeatureStorage(this).mat_node_id_;
  auto* result = mat.get_node(id);
  Assert(result != nullptr);
  return *result;
}

template <typename CRTP, typename Tag>
MAT::Node& FeatureMutableView<MATConversion, CRTP, Tag>::GetMATNode() const {
  auto& node = static_cast<const CRTP&>(*this);
  auto& mat = node.GetDAG().GetMAT();
  size_t id = GetFeatureStorage(this).mat_node_id_;
  auto* result = mat.get_node(id);
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
MAT::Tree& ExtraFeatureMutableView<MATConversion, CRTP>::GetMAT() const {
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

template <typename DAG>
static void mat_from_dag_helper(typename DAG::NodeView dag_node,
                                MAT::Node* mat_par_node, MAT::Tree& new_tree) {
  mat_par_node->children.reserve(dag_node.GetCladesCount());
  for (auto clade : dag_node.GetClades()) {
    Assert(clade.size() == 1);
    typename DAG::EdgeView edge = *clade.begin();
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
    mat_from_dag_helper<DAG>(edge.GetChild(), node, new_tree);
  }
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
  mat_from_dag_helper<std::decay_t<decltype(dag)>>(root_node, mat_root_node, tree);
  return tree;
}
