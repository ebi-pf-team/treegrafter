#!/usr/bin/env bash

set -e

cd "$(dirname "$0")"

if [ "$1" == "prepare" ]; then
  python3 treegrafter.py prepare "${@:2}"
elif [ "$1" == "search" ]; then
  fasta="$2"
  datadir="$3"
  hmmsearchout=$(mktemp)
  hmmsearch "${datadir}"/famhmm/binHmm "${fasta}" > "${hmmsearchout}"
  python3 treegrafter.py run "${fasta}" "${hmmsearchout}" "${datadir}" "${@:4}"
fi
