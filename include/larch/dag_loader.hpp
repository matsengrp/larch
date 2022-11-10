/**
 * Functions for loading and writing mutation annotated dags.
 */
#pragma once

#include <fstream>
#include <vector>
#include <unordered_map>
#include <string_view>

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

#include "larch/merge/merge.hpp"

[[nodiscard]] inline MADAGStorage LoadDAGFromProtobuf(std::string_view path);

[[nodiscard]] inline MADAGStorage LoadTreeFromProtobuf(
    std::string_view path, std::string_view reference_sequence);

[[nodiscard]] inline MADAGStorage LoadDAGFromJson(std::string_view path);

[[nodiscard]] inline std::string LoadReferenceSequence(std::string_view path);

template <typename DAG>
void StoreDAGToProtobuf(DAG dag, std::string_view path);

void inline StoreTreeToProtobuf(MADAG dag, std::string_view path);

void inline MADAGToDOT(MADAG dag, std::ostream& out);

#include "larch/impl/dag_loader_impl.hpp"