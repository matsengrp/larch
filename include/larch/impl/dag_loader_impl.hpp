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