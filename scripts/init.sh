#!/usr/bin/env bash
# Initialize larch-build.env with build settings.
#
# Settings come from environment variables, falling back to
# existing larch-build.env values, then defaults.
#
# Usage:
#   pixi run init
#   LARCH_BUILD_DIR=build-dev LARCH_USE_USHER=OFF pixi run init
#
# After init, edit larch-build.env as needed, then run: pixi run configure

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

# Initialize git submodules
echo "Updating git submodules..."
git -C "$PROJECT_DIR" submodule update --init --recursive
CONFIG_FILE="$PROJECT_DIR/larch-build.env"

# Source existing config as defaults (if present)
if [[ -f "$CONFIG_FILE" ]]; then
    source "$CONFIG_FILE"
fi

# Apply env vars over existing config, fall back to defaults
LARCH_BUILD_DIR="${LARCH_BUILD_DIR:-build}"
LARCH_BUILD_TYPE="${LARCH_BUILD_TYPE:-RelWithDebInfo}"
LARCH_USE_USHER="${LARCH_USE_USHER:-ON}"
LARCH_USE_NETAM="${LARCH_USE_NETAM:-OFF}"
LARCH_LIBTORCH_PATH="${LARCH_LIBTORCH_PATH:-auto}"
LARCH_CMAKE_EXTRA="${LARCH_CMAKE_EXTRA:-}"
LARCH_BUILD_JOBS="${LARCH_BUILD_JOBS:-$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)}"

# Write config file
cat > "$CONFIG_FILE" <<EOF
# Larch build configuration
# Edit this file, then run: pixi run configure

# Build directory name
LARCH_BUILD_DIR=$LARCH_BUILD_DIR

# CMake build type: Debug, Release, RelWithDebInfo, MinSizeRel
LARCH_BUILD_TYPE=$LARCH_BUILD_TYPE

# Enable matOptimize integration (ON/OFF)
LARCH_USE_USHER=$LARCH_USE_USHER

# Enable ML scoring models via libtorch (ON/OFF)
LARCH_USE_NETAM=$LARCH_USE_NETAM

# Path to libtorch (used when LARCH_USE_NETAM=ON)
#   "auto" = use pixi/conda environment (default)
#   custom path = e.g. /home/user/libtorch or deps/libtorch
LARCH_LIBTORCH_PATH=$LARCH_LIBTORCH_PATH

# Additional cmake flags. Supported options:
#   -DUSE_MAT_VIEW=OFF        Use MAT conversion instead of MAT view
#   -DUSE_ASAN=yes            Enable address/UB sanitizer
#   -DUSE_TSAN=yes            Enable thread sanitizer
#   -DDISABLE_PARALLELISM=yes Disable all parallelism (debugging)
#   -DKEEP_ASSERTS=ON         Keep asserts in release builds
#   -DUSE_CPPTRACE=ON         Readable backtraces on exceptions
#   -DUSE_SYSTEM_TBB=yes      Use system TBB instead of fetching from source
# Example: LARCH_CMAKE_EXTRA=-DUSE_ASAN=yes -DKEEP_ASSERTS=ON
LARCH_CMAKE_EXTRA=$LARCH_CMAKE_EXTRA

# Number of parallel build jobs
LARCH_BUILD_JOBS=$LARCH_BUILD_JOBS
EOF

echo "=== Larch build configuration ==="
cat "$CONFIG_FILE"
echo "================================="
echo "Config written to $CONFIG_FILE"
echo "Edit as needed, then run: pixi run configure"
