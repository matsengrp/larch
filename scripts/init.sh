#!/usr/bin/env bash
# Initialize larch-build.env from template.
# Environment variables override template defaults.
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
TEMPLATE="$SCRIPT_DIR/larch-build.env.template"

if [[ -f "$CONFIG_FILE" ]]; then
    echo "larch-build.env already exists. Remove it to re-initialize from template."
    echo "Current config:"
    cat "$CONFIG_FILE"
else
    # Copy template, then apply any environment variable overrides
    cp "$TEMPLATE" "$CONFIG_FILE"
    for var in LARCH_BUILD_DIR LARCH_BUILD_TYPE LARCH_USE_USHER LARCH_USE_NETAM \
               LARCH_LIBTORCH_PATH LARCH_PROTOBUF_PATH LARCH_CMAKE_EXTRA LARCH_BUILD_JOBS; do
        if [[ -n "${!var:-}" ]]; then
            sed -i "s|^${var}=.*|${var}=${!var}|" "$CONFIG_FILE"
            echo "Override: ${var}=${!var}"
        fi
    done
    echo "Config created at $CONFIG_FILE"
    cat "$CONFIG_FILE"
fi
