#pragma once

#include <string_view>
#include <unordered_map>

#include "merge.hpp"

[[nodiscard]] HistoryDAG LoadHistoryDAGFromProtobufGZ(
    std::string_view path, std::string& ref_seq, std::vector<Mutations>& mutations);

[[nodiscard]] HistoryDAG LoadTreeFromProtobufGZ(std::string_view path,
                                                std::vector<Mutations>& mutations);

[[nodiscard]] std::string LoadRefseqFromJsonGZ(std::string_view path);

[[nodiscard]] HistoryDAG LoadHistoryDAGFromJsonGZ(std::string_view path,
                                                  std::string& refseq);

[[nodiscard]] std::vector<NodeLabel> LoadLabelsJsonGZ(std::string_view path);

void StoreDAGToProtobuf(const HistoryDAG& dag, std::string_view reference_sequence,
                        const ConcurrentUnorderedMap<NodeLabel, NodeId>& labels,
                        std::string_view path);
