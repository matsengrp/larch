#pragma once

#include <vector>
#include <tuple>
#include <map>
#include <type_traits>

struct NodeId;
struct EdgeId;
struct CladeIdx;

template <typename... Features>
class DefaultNodeStorage;

template <typename Storage>
class DefaultNodesContainer;

template <typename... Features>
class DefaultEdgeStorage;

template <typename Storage>
class DefaultEdgesContainer;

template <typename NodesContainer, typename EdgesContainer, typename... Features>
class DefaultDAGStorage;

template <typename Storage, typename... Features>
class DAGView;
template <typename DAG, typename... Features>
class NodeView;
template <typename DAG, typename... Features>
class EdgeView;

template <typename Feature, typename View>
class FeatureReader;
template <typename Feature, typename View>
class FeatureWriter;

#define DAG_DECLARATIONS

#include "larch/dag/dag_common.hpp"
#include "larch/dag/node_storage.hpp"
#include "larch/dag/nodes_container.hpp"
#include "larch/dag/edge_storage.hpp"
#include "larch/dag/edges_container.hpp"
#include "larch/dag/dag_storage.hpp"
#include "larch/dag/dag_view.hpp"
#include "larch/dag/node_view.hpp"
#include "larch/dag/edge_view.hpp"

#undef DAG_DECLARATIONS

#define DAG_DEFINITIONS

#include "larch/impl/dag/dag_common_impl.hpp"
#include "larch/impl/dag/node_storage_impl.hpp"
#include "larch/impl/dag/nodes_container_impl.hpp"
#include "larch/impl/dag/edge_storage_impl.hpp"
#include "larch/impl/dag/edges_container_impl.hpp"
#include "larch/impl/dag/dag_storage_impl.hpp"
#include "larch/impl/dag/dag_view_impl.hpp"
#include "larch/impl/dag/node_view_impl.hpp"
#include "larch/impl/dag/edge_view_impl.hpp"

#undef DAG_DEFINITIONS
