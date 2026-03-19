#!/usr/bin/env bash
# Run larch-test with the build config. Passes all arguments to larch-test.
#
# Usage:
#   bash scripts/run-tests.sh -tag slow    # exclude slow tests
#   bash scripts/run-tests.sh              # run all tests

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/load-config.sh"

mkdir -p "$LARCH_BUILD_DIR/tests"
ln -sfn ../../data "$LARCH_BUILD_DIR/tests/data"
cd "$LARCH_BUILD_DIR/tests"
exec ../bin/larch-test "$@"
