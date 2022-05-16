#pragma once

#include <vector>
#include <functional>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <tuple>
#include <cassert>

#include "history_dag.hpp"

using CompactGenome = std::map<size_t, char>;
using LeafSet = std::set<std::set<CompactGenome>>;
using NodeLabel = std::pair<CompactGenome, LeafSet>;
using EdgeLabel = std::tuple<NodeLabel, NodeLabel, CladeIdx>;

class Merge;

class Merge {
public:
    inline Merge(const std::string& refseq,
        std::vector<std::reference_wrapper<const HistoryDAG>>&& trees,
        const std::vector<std::vector<NodeLabel>>& labels);

    inline void Run();

    inline HistoryDAG& GetResult();
    inline const HistoryDAG& GetResult() const;

    inline const NodeLabel& GetNodeLabel(size_t tree_idx, NodeId node_id);
    inline NodeId GetResultNode(const NodeLabel& label);
    inline void MakeResultEdge(size_t tree_idx, EdgeId edge_id);
    

private:
    const std::string refseq_; 
    std::vector<std::reference_wrapper<const HistoryDAG>> trees_;
    const std::vector<std::vector<NodeLabel>>& labels_;
    std::map<NodeLabel, NodeId> result_nodes_;
    std::set<EdgeLabel> result_edges_;
    HistoryDAG result_;
};

#include "impl/merge_impl.hpp"
