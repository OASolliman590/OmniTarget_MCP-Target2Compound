#!/usr/bin/env bash
set -euo pipefail

echo "Preparing Ouroboros model directory in third_party/Ouroboros/models ..."

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
REPO_DIR="$ROOT_DIR/third_party/Ouroboros"
DEST_DIR="$REPO_DIR/models"
mkdir -p "$DEST_DIR"

echo "This script does not auto-download Ouroboros models."
echo "Per upstream README, download models from ZhangLab WebPage:"
echo "  https://zhanglab.comp.nus.edu.sg/Ouroboros/"
echo "Then place the chosen model directory (e.g., Ouroboros_M1c) under:"
echo "  $DEST_DIR/"
echo
echo "Example final structure:"
cat <<'TREE'
third_party/Ouroboros/
  models/
    Ouroboros_M1c/
      config.json
      encoder_weights.pth
      decoder_weights.pth
      feat_stat.csv
      ...
TREE

echo "Creating placeholder directory: $DEST_DIR/Ouroboros_M1c"
mkdir -p "$DEST_DIR/Ouroboros_M1c"
touch "$DEST_DIR/Ouroboros_M1c/.place_model_files_here"
echo "✅ Placeholder created. Copy the actual model files into this folder."

