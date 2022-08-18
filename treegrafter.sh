#!/usr/bin/env bash

set -e

cd "$(dirname "$0")"

if [ "$1" == "prepare" ]; then
  python3 treegrafter.py prepare "${@:2}"
elif [ "$1" == "search" ]; then
  fasta="$2"
  datadir="$3"
  hmmsearchout=$(mktemp)
  treegrafterout="$4"
  hmmsearch -Z 65000000 -E 0.001 --domE 0.00000001 --incdomE 0.00000001 --notextw "$datadir"/famhmm/binHmm "$fasta" > "$hmmsearchout"
  python3 treegrafter.py run -e 0.00000001 -o "$treegrafterout" "$fasta" "$hmmsearchout" "$datadir"
fi
