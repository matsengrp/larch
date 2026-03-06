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

# Resolve protobuf path
PROTOBUF_ARGS=""
if [[ "$LARCH_PROTOBUF_PATH" == "auto" ]]; then
  # Ensure cmake finds protobuf headers from pixi/conda, not the system
  if [[ -n "${CONDA_PREFIX:-}" && -d "$CONDA_PREFIX/include/google/protobuf" ]]; then
    PROTOBUF_ARGS="-DProtobuf_INCLUDE_DIR=$CONDA_PREFIX/include"
    echo "Using protobuf from conda environment: $CONDA_PREFIX"
  fi
else
  # Resolve relative paths against project dir
  if [[ "$LARCH_PROTOBUF_PATH" != /* ]]; then
    PROTOBUF_ABS="$PROJECT_DIR/$LARCH_PROTOBUF_PATH"
  else
    PROTOBUF_ABS="$LARCH_PROTOBUF_PATH"
  fi
  if [[ ! -d "$PROTOBUF_ABS/include/google/protobuf" ]]; then
    echo "Error: protobuf headers not found at $PROTOBUF_ABS/include/google/protobuf" >&2
    echo "Set LARCH_PROTOBUF_PATH to a prefix containing include/ and lib/ directories." >&2
    exit 1
  fi
  PROTOBUF_ARGS="-DProtobuf_INCLUDE_DIR=$PROTOBUF_ABS/include -DProtobuf_LIBRARY=$PROTOBUF_ABS/lib/libprotobuf.so -DProtobuf_PROTOC_EXECUTABLE=$PROTOBUF_ABS/bin/protoc"
  echo "Using protobuf from: $PROTOBUF_ABS"
fi

# Save a copy of the build config into the build directory for reference
mkdir -p "$PROJECT_DIR/$LARCH_BUILD_DIR"
cp "$SCRIPT_DIR/../larch-build.env" "$PROJECT_DIR/$LARCH_BUILD_DIR/larch-build.env"

cmake -S "$PROJECT_DIR" -B "$PROJECT_DIR/$LARCH_BUILD_DIR" \
  -DCMAKE_BUILD_TYPE="$LARCH_BUILD_TYPE" \
  -DUSE_USHER="$LARCH_USE_USHER" \
  -DUSE_NETAM="$LARCH_USE_NETAM" \
  -DTBB_DISABLE_HWLOC_AUTOMATIC_SEARCH=ON \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
  $PROTOBUF_ARGS \
  $CMAKE_PREFIX_PATH_ARG \
  $LARCH_CMAKE_EXTRA
