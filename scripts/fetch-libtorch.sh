#!/usr/bin/env bash
# Download and extract the libtorch C++ binary distribution.
# Installs to deps/libtorch/ by default.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

LIBTORCH_VERSION="${LIBTORCH_VERSION:-2.10.0}"
LIBTORCH_DIR="$PROJECT_DIR/deps/libtorch"

if [[ -d "$LIBTORCH_DIR" ]]; then
    echo "libtorch already exists at $LIBTORCH_DIR"
    echo "Remove it to re-download: rm -rf $LIBTORCH_DIR"
    exit 0
fi

# Detect platform
OS="$(uname -s)"
ARCH="$(uname -m)"

case "$OS" in
    Linux)
        URL="https://download.pytorch.org/libtorch/cpu/libtorch-shared-with-deps-${LIBTORCH_VERSION}%2Bcpu.zip"
        ;;
    Darwin)
        if [[ "$ARCH" == "arm64" ]]; then
            URL="https://download.pytorch.org/libtorch/cpu/libtorch-macos-arm64-${LIBTORCH_VERSION}.zip"
        else
            URL="https://download.pytorch.org/libtorch/cpu/libtorch-macos-x86_64-${LIBTORCH_VERSION}.zip"
        fi
        ;;
    *)
        echo "Error: Unsupported platform: $OS" >&2
        exit 1
        ;;
esac

echo "Downloading libtorch ${LIBTORCH_VERSION} for ${OS}/${ARCH}..."
echo "URL: $URL"

TMPFILE="$(mktemp /tmp/libtorch-XXXXXX.zip)"
trap "rm -f '$TMPFILE'" EXIT

curl -L -o "$TMPFILE" "$URL"

echo "Extracting to $LIBTORCH_DIR..."
mkdir -p "$PROJECT_DIR/deps"
unzip -q "$TMPFILE" -d "$PROJECT_DIR/deps"

echo "libtorch installed to $LIBTORCH_DIR"
