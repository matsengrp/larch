#pragma once

#include <string_view>
#include <unordered_map>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wparentheses"
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic ignored "-Wclass-memaccess"
#pragma GCC diagnostic ignored "-Wdeprecated-copy"
#pragma GCC diagnostic ignored "-Wformat"
#pragma GCC diagnostic ignored "-Wextra"
#pragma GCC diagnostic ignored "-Wpedantic"
#include "src/matOptimize/mutation_annotated_tree.hpp"
#pragma GCC diagnostic pop

#include "merge.hpp"

[[nodiscard]] MADAG LoadDAGFromProtobuf(std::string_view path);

[[nodiscard]] MADAG LoadTreeFromProtobuf(std::string_view path);

[[nodiscard]] MADAG LoadDAGFromJson(std::string_view path);

void StoreDAGToProtobuf(const DAG& dag, std::string_view reference_sequence,
                        const std::vector<EdgeMutations>& edge_parent_mutations,
                        std::string_view path);

void MADAGToDOT(const MADAG& dag, std::ostream& out);

[[nodiscard]] MADAG LoadDAGFromUsher(const Mutation_Annotated_Tree::Tree* tree,
                                     std::string_view reference_sequence);

Mutation_Annotated_Tree::Tree StoreDAGToUsher(const MADAG& dag);
