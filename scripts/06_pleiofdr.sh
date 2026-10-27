#!/usr/bin/env bash
set -euo pipefail
# Usage: 06_pleiofdr.sh /path/to/config.txt /path/to/patched/run_pair.m
config=${1:?config file}; runner=${2:?patched pleioFDR runner}
: "${MATLAB_BIN:?}" "${PLEIOFDR_ROOT:?}"
"$MATLAB_BIN" -batch "addpath('$PLEIOFDR_ROOT'); config='$config'; run('$runner')"
