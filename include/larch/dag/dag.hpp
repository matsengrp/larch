#pragma once

#include <vector>
#include <tuple>
#include <type_traits>

#include <tbb/concurrent_unordered_set.h>

#include "larch/common.hpp"

#define DAG_DECLARATIONS
#include "larch/dag/dag_common.hpp"
#include "larch/dag/element_storage.hpp"
#include "larch/dag/elements_container.hpp"
#include "larch/dag/dag_storage.hpp"
#include "larch/dag/element_view.hpp"
#include "larch/dag/dag_view.hpp"
#include "larch/dag/extend.hpp"
#include "larch/dag/deduplicate.hpp"
#include "larch/dag/neighbors.hpp"
#include "larch/dag/endpoints.hpp"
#include "larch/dag/connections.hpp"
#undef DAG_DECLARATIONS

#define DAG_DEFINITIONS
#include "larch/impl/dag/dag_common_impl.hpp"
#include "larch/impl/dag/element_storage_impl.hpp"
#include "larch/impl/dag/elements_container_impl.hpp"
#include "larch/impl/dag/dag_storage_impl.hpp"
#include "larch/impl/dag/element_view_impl.hpp"
#include "larch/impl/dag/dag_view_impl.hpp"
#include "larch/impl/dag/extend_impl.hpp"
#include "larch/impl/dag/deduplicate_impl.hpp"
#include "larch/impl/dag/neighbors_impl.hpp"
#include "larch/impl/dag/endpoints_impl.hpp"
#include "larch/impl/dag/connections_impl.hpp"
#undef DAG_DEFINITIONS
