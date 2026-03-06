#!/usr/bin/env bash
# Source larch-build.env and export settings.
# This script is meant to be sourced (not executed) by other tasks.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="$SCRIPT_DIR/../larch-build.env"

if [[ -f "$CONFIG_FILE" ]]; then
    echo "Loading config from $CONFIG_FILE"
    source "$CONFIG_FILE"
else
    echo "Error: No larch-build.env found. Run 'pixi run init' first." >&2
    exit 1
fi

# Default LARCH_BUILD_JOBS to number of CPU cores if empty
export LARCH_BUILD_JOBS="${LARCH_BUILD_JOBS:-$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)}"

export LARCH_BUILD_DIR LARCH_BUILD_TYPE LARCH_USE_USHER LARCH_USE_NETAM
export LARCH_LIBTORCH_PATH LARCH_PROTOBUF_PATH LARCH_CMAKE_EXTRA

echo "=== Build settings ==="
echo "  LARCH_BUILD_DIR=$LARCH_BUILD_DIR"
echo "  LARCH_BUILD_TYPE=$LARCH_BUILD_TYPE"
echo "  LARCH_USE_USHER=$LARCH_USE_USHER"
echo "  LARCH_USE_NETAM=$LARCH_USE_NETAM"
echo "  LARCH_LIBTORCH_PATH=$LARCH_LIBTORCH_PATH"
echo "  LARCH_PROTOBUF_PATH=$LARCH_PROTOBUF_PATH"
echo "  LARCH_BUILD_JOBS=$LARCH_BUILD_JOBS"
echo "======================"
