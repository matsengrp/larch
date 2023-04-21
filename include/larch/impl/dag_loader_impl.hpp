#include <fstream>
#include <vector>
#include <unordered_map>

#include <sys/stat.h>
#include <fcntl.h>

#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/gzip_stream.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstack-usage="
#include "nlohmann/json.hpp"
#pragma GCC diagnostic pop

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/copy.hpp>

#include "larch/dag_loader.hpp"
#include "dag.pb.h"
#include "parsimony.pb.h"
#include "larch/newick.hpp"

inline int32_t EncodeBasePB(char base) {
  switch (base) {
    case 'A':
      return 0;
    case 'C':
      return 1;
    case 'G':
      return 2;
    case 'T':
      return 3;
    default:
      Fail("Invalid base");
  };
}

template <typename DAG>
void StoreDAGToProtobuf(DAG dag, std::string_view path) {
  dag.AssertUA();
  ProtoDAG::data data;

  data.set_reference_seq(dag.GetReferenceSequence());

  for (size_t i = 0; i < dag.GetNodesCount(); ++i) {
    auto* proto_node = data.add_node_names();
    proto_node->set_node_id(static_cast<int64_t>(i));
  }

  for (typename DAG::EdgeView edge : dag.GetEdges()) {
    auto* proto_edge = data.add_edges();
    proto_edge->set_edge_id(static_cast<int64_t>(edge.GetId().value));
    proto_edge->set_parent_node(static_cast<int64_t>(edge.GetParentId().value));
    proto_edge->set_child_node(static_cast<int64_t>(edge.GetChildId().value));
    proto_edge->set_parent_clade(static_cast<int64_t>(edge.GetClade().value));
    for (auto [pos, nucs] : edge.GetEdgeMutations()) {
      auto* proto_mut = proto_edge->add_edge_mutations();
      proto_mut->set_position(static_cast<int32_t>(pos.value));
      proto_mut->set_par_nuc(EncodeBasePB(nucs.first));
      proto_mut->add_mut_nuc(EncodeBasePB(nucs.second));
    }
  }

  std::ofstream file{std::string{path}};
  data.SerializeToOstream(&file);
}

template <typename DAG>
void StoreTreeToProtobuf(DAG dag, std::string_view path) {
  dag.AssertUA();
  Assert(dag.IsTree());

  Parsimony::data data;
  using Node = typename DAG::NodeView;
  std::string newick;
  auto to_newick = [&](auto& self, Node node) -> void {
    if (not node.IsLeaf()) {
      newick += '(';
    }
    size_t clade_idx = 0;
    for (auto clade : node.GetClades()) {
      Assert(clade.size() == 1);
      Node i = (*clade.begin()).GetChild();
      if (i.IsLeaf()) {
        if (i.GetSampleId()) {
          newick += *i.GetSampleId();
        } else {
          newick += "unknown_leaf_";
          newick += std::to_string(i.GetId().value);
        }
      } else {
        self(self, i);
        newick += "inner_";
        newick += std::to_string(i.GetId().value);
      }
      if (++clade_idx < node.GetCladesCount()) {
        newick += ',';
      }
    }
    if (not node.IsLeaf()) {
      newick += ')';
    }
  };
  to_newick(to_newick, dag.GetRoot().GetFirstChild().GetChild());
  newick += ';';
  data.set_newick(newick);
  using Edge = typename DAG::EdgeView;
  auto store_mutations = [](auto& self, Edge edge, Parsimony::data& result,
                            std::string_view ref_seq) -> void {
    auto* proto = result.add_node_mutations();
    for (auto [pos, mut] : edge.GetEdgeMutations()) {
      auto* proto_mut = proto->add_mutation();
      proto_mut->set_position(static_cast<int32_t>(pos.value));
      proto_mut->set_ref_nuc(EncodeBasePB(ref_seq.at(pos.value - 1)));
      proto_mut->set_par_nuc(EncodeBasePB(mut.first));
      proto_mut->add_mut_nuc(EncodeBasePB(mut.second));
      proto_mut->set_chromosome("leaf_0");
    }
    for (Edge child : edge.GetChild().GetChildren()) {
      self(self, child, result, ref_seq);
    }
  };
  store_mutations(store_mutations, dag.GetRoot().GetFirstChild(), data,
                  dag.GetReferenceSequence());

  std::ofstream file{std::string{path}};
  data.SerializeToOstream(&file);
}

template <typename Edge>
static std::string EdgeMutationsToString(Edge edge) {
  std::string result;
  size_t count = 0;
  for (auto [pos, muts] : edge.GetEdgeMutations()) {
    result += muts.first;
    result += std::to_string(pos.value);
    result += muts.second;
    result += ++count % 3 == 0 ? "\\n" : " ";
  }
  return result;
}

template <typename Node>
static std::string CompactGenomeToString(Node node) {
  if (node.IsUA()) {
    return std::to_string(node.GetId().value);
  }
  std::string result = std::to_string(node.GetId().value);
  if (node.HasChangedTopology()) {
    result += " [X]";
  }
  result += "\\n";
  size_t count = 0;
  for (auto [pos, base] : node.Const().GetCompactGenome()) {
    result += std::to_string(pos.value);
    result += base;
    result += ++count % 3 == 0 ? "\\n" : " ";
  }
  return result;
}

template <typename DAG>
void MADAGToDOT(DAG dag, std::ostream& out) {
  out << "digraph {\n";
  out << "  forcelabels=true\n";
  out << "  nodesep=1.0\n";
  out << "  ranksep=2.0\n";
  out << "  ratio=1.0\n";
  for (auto edge : dag.Const().GetEdges()) {
    out << "  \"" << CompactGenomeToString(edge.GetParent()) << "\" -> \""
        << CompactGenomeToString(edge.GetChild()) << "\"";
    out << "[ xlabel=\"";
    out << EdgeMutationsToString(edge);
    out << "\" ]";
    out << "\n";
  }
  out << "}\n";
}

template <typename DAG>
void FragmentToDOT(DAG dag, const std::vector<EdgeId>& edges, std::ostream& out) {
  out << "digraph {\n";
  out << "  forcelabels=true\n";
  out << "  nodesep=1.0\n";
  out << "  ranksep=2.0\n";
  out << "  ratio=1.0\n";
  for (auto edge : edges | Transform::ToEdges(dag)) {
    out << "  \"" << CompactGenomeToString(edge.GetParent()) << "\" -> \""
        << CompactGenomeToString(edge.GetChild()) << "\"";
    out << "[ xlabel=\"";
    out << EdgeMutationsToString(edge);
    out << "\" ]";
    out << "\n";
  }
  out << "}\n";
}
