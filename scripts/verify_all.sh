#!/usr/bin/env bash
# Rebuild VERIFICATION.txt: the test suite and every numbered script, run in
# order with the project's own interpreter, so the log records the environment
# the thesis numbers come from. Also regenerates every figure.
#
#   bash scripts/verify_all.sh            (about ten minutes; script 06 is slow)
#
# Uses .venv/bin/python if it exists, else whatever `python` is. The versions it
# ran with are printed at the top of the log; compare them with
# requirements-lock.txt before trusting a residual that differs from the thesis.
set -euo pipefail
cd "$(dirname "$0")/.."
PY=python
[ -x .venv/bin/python ] && PY=.venv/bin/python
OUT=VERIFICATION.txt
{
  echo "VERIFICATION LOG"
  # "-dirty" means the run used uncommitted changes on top of that commit.
  echo "schwinger-n4  |  generated $(date -u '+%Y-%m-%d %H:%M UTC')  |  commit $(git describe --always --dirty 2>/dev/null || echo none)"
  $PY -c "import sys, numpy, scipy, qiskit, qiskit_aer, matplotlib; print(f'python {sys.version.split()[0]}  numpy {numpy.__version__}  scipy {scipy.__version__}  qiskit {qiskit.__version__}  qiskit-aer {qiskit_aer.__version__}  matplotlib {matplotlib.__version__}')"
  echo
  echo "############ TEST SUITE ############"
  $PY -m pytest tests/ -q 2>&1 | tail -3
  for s in scripts/0*.py; do
    echo
    echo "############ $s ############"
    $PY "$s" 2>&1 | sed "s|$(pwd)/||g"
  done
} > "$OUT"
echo "wrote $OUT"
