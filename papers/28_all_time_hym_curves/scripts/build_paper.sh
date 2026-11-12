#!/usr/bin/env bash
set -euo pipefail
paper28_dir="$(cd "$(dirname "$0")/.." && pwd)"
latexmk -pdf -interaction=nonstopmode -halt-on-error -cd "$paper28_dir/main.tex"
