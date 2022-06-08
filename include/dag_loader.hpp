#pragma once

#include <string_view>
#include <unordered_map>

#include "merge.hpp"

[[nodiscard]] MADAG LoadDAGFromProtobuf(std::string_view path);

[[nodiscard]] MADAG LoadTreeFromProtobuf(std::string_view path);

[[nodiscard]] MADAG LoadDAGFromJson(std::string_view path);

[[nodiscard]] std::vector<CompactGenome> LoadCompactGenomesJson(std::string_view path);

void StoreDAGToProtobuf(const DAG& dag, std::string_view reference_sequence,
                        const std::vector<EdgeMutations>& edge_parent_mutations,
                        std::string_view path);

void EdgeMutationsDAGToDOT(const MADAG& dag, std::ostream& out);

void CompactGenomeDAGToDOT(const MADAG& dag,
                           const std::vector<CompactGenome>& compact_genomes,
                           std::ostream& out);
