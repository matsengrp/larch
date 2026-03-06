#!/usr/bin/env bash
# Initialize larch-build.env from template.
#
# Usage:
#   pixi run init
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
    cp "$TEMPLATE" "$CONFIG_FILE"
    echo "Config created at $CONFIG_FILE"
    cat "$CONFIG_FILE"
fi
