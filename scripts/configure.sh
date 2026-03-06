#!/usr/bin/env bash
# Run cmake using settings from larch-build.env.
# Run 'pixi run init' first to create the config file.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

source "$SCRIPT_DIR/load-config.sh"

# Build CMAKE_PREFIX_PATH for libtorch when USE_NETAM is enabled
CMAKE_PREFIX_PATH_ARG=""
if [[ "$LARCH_USE_NETAM" == "ON" || "$LARCH_USE_NETAM" == "yes" ]]; then
  if [[ "$LARCH_LIBTORCH_PATH" == "auto" ]]; then
    if [[ -n "${CONDA_PREFIX:-}" ]]; then
      # Find torch cmake config (may be in site-packages or directly in prefix)
      TORCH_CMAKE_DIR=$(find "$CONDA_PREFIX" -name "TorchConfig.cmake" -printf '%h\n' 2>/dev/null | head -1)
      if [[ -z "$TORCH_CMAKE_DIR" ]]; then
        echo "Error: pytorch-cpu is installed but TorchConfig.cmake not found in $CONDA_PREFIX" >&2
        exit 1
      fi
      TORCH_PREFIX=$(dirname $(dirname $(dirname "$TORCH_CMAKE_DIR")))
      CMAKE_PREFIX_PATH_ARG="-DCMAKE_PREFIX_PATH=$CONDA_PREFIX;$TORCH_PREFIX"
      echo "Using libtorch from conda environment: $TORCH_PREFIX"
    else
      echo "Error: LARCH_LIBTORCH_PATH is 'auto' but CONDA_PREFIX is not set." >&2
      echo "Set LARCH_LIBTORCH_PATH in larch-build.env to your libtorch install path." >&2
      exit 1
    fi
  else
    # Resolve relative paths against project dir
    if [[ "$LARCH_LIBTORCH_PATH" != /* ]]; then
      LIBTORCH_ABS="$PROJECT_DIR/$LARCH_LIBTORCH_PATH"
    else
      LIBTORCH_ABS="$LARCH_LIBTORCH_PATH"
    fi
    if [[ ! -d "$LIBTORCH_ABS" ]]; then
      echo "Error: libtorch not found at $LIBTORCH_ABS" >&2
      echo "Run 'pixi run fetch-libtorch' or set LARCH_LIBTORCH_PATH in larch-build.env." >&2
      exit 1
    fi
    CMAKE_PREFIX_PATH_ARG="-DCMAKE_PREFIX_PATH=$LIBTORCH_ABS"
    echo "Using libtorch from: $LIBTORCH_ABS"
  fi
fi

# Ensure cmake finds protobuf headers from pixi/conda, not the system
PROTOBUF_INCLUDE_ARG=""
if [[ -n "${CONDA_PREFIX:-}" && -d "$CONDA_PREFIX/include/google/protobuf" ]]; then
  PROTOBUF_INCLUDE_ARG="-DProtobuf_INCLUDE_DIR=$CONDA_PREFIX/include"
fi

cmake -S "$PROJECT_DIR" -B "$PROJECT_DIR/$LARCH_BUILD_DIR" \
  -DCMAKE_BUILD_TYPE="$LARCH_BUILD_TYPE" \
  -DUSE_USHER="$LARCH_USE_USHER" \
  -DUSE_NETAM="$LARCH_USE_NETAM" \
  -DTBB_DISABLE_HWLOC_AUTOMATIC_SEARCH=ON \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
  $PROTOBUF_INCLUDE_ARG \
  $CMAKE_PREFIX_PATH_ARG \
  $LARCH_CMAKE_EXTRA
