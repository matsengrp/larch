#!/usr/bin/env bash
# Source larch-build.env if it exists, otherwise use defaults.
# This script is meant to be sourced (not executed) by other tasks.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="$SCRIPT_DIR/../larch-build.env"

if [[ -f "$CONFIG_FILE" ]]; then
    echo "Loading config from $CONFIG_FILE"
    source "$CONFIG_FILE"
else
    echo "Error: No larch-build.env found. Run 'pixi run init' first." >&2
    exit 1
fi

export LARCH_BUILD_DIR="${LARCH_BUILD_DIR:-build}"
export LARCH_BUILD_TYPE="${LARCH_BUILD_TYPE:-Debug}"
export LARCH_USE_USHER="${LARCH_USE_USHER:-ON}"
export LARCH_USE_NETAM="${LARCH_USE_NETAM:-OFF}"
export LARCH_LIBTORCH_PATH="${LARCH_LIBTORCH_PATH:-deps/libtorch}"
export LARCH_CMAKE_EXTRA="${LARCH_CMAKE_EXTRA:-}"
export LARCH_BUILD_JOBS="${LARCH_BUILD_JOBS:-$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)}"
