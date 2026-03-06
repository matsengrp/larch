#!/usr/bin/env bash
# Run cmake using settings from larch-build.env.
# Run 'pixi run init' first to create the config file.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

source "$SCRIPT_DIR/load-config.sh"

echo "=== Configuring with ==="
cat "$PROJECT_DIR/larch-build.env"
echo "========================"
echo ""

# Build CMAKE_PREFIX_PATH for libtorch when USE_NETAM is enabled
CMAKE_PREFIX_PATH_ARG=""
if [[ "$LARCH_USE_NETAM" == "ON" || "$LARCH_USE_NETAM" == "yes" ]]; then
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

cmake -S "$PROJECT_DIR" -B "$PROJECT_DIR/$LARCH_BUILD_DIR" \
  -DCMAKE_BUILD_TYPE="$LARCH_BUILD_TYPE" \
  -DUSE_USHER="$LARCH_USE_USHER" \
  -DUSE_NETAM="$LARCH_USE_NETAM" \
  -DTBB_DISABLE_HWLOC_AUTOMATIC_SEARCH=ON \
  $CMAKE_PREFIX_PATH_ARG \
  $LARCH_CMAKE_EXTRA
