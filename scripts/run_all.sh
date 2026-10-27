#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
for input in data/PA_7071.xlsx data/PA_7071_BioProject.xlsx; do
  if [[ ! -f "$input" ]]; then
    echo "Missing required input: $input" >&2
    echo "See README.md: input data are not included in this code-only release." >&2
    exit 1
  fi
done
mkdir -p results
python scripts/run_core_analyses.py --input data/PA_7071.xlsx --output results
python scripts/run_country_sensitivity.py --input data/PA_7071.xlsx --output results
python scripts/run_bioproject_robustness.py
