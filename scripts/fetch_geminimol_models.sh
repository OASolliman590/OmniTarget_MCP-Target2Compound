#!/usr/bin/env bash
set -euo pipefail

echo "Fetching GeminiMol model files into third_party/GeminiMol/models/GeminiMol ..."

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
DEST_DIR="$ROOT_DIR/third_party/GeminiMol/models/GeminiMol"
mkdir -p "$DEST_DIR"

cd "$DEST_DIR"

echo "Destination: $(pwd)"

have_gdown=0
have_hf=0
if command -v gdown >/dev/null 2>&1; then have_gdown=1; fi
if command -v huggingface-cli >/dev/null 2>&1; then have_hf=1; fi

if [ $have_gdown -eq 1 ] || python -c "import gdown" >/dev/null 2>&1; then
  echo "Using gdown to fetch GeminiMol.zip (official link from README)"
  echo "If this fails, try manual download via browser."
  if command -v gdown >/dev/null 2>&1; then
    gdown https://drive.google.com/uc?id=1IgpkabylSJ0aIwUjdIWOC5ID51pcucw5 -O GeminiMol.zip || {
      echo "gdown failed. You may need to authenticate or fetch manually."; exit 1;
    }
  else
    python -m gdown https://drive.google.com/uc?id=1IgpkabylSJ0aIwUjdIWOC5ID51pcucw5 -O GeminiMol.zip || {
      echo "python -m gdown failed. You may need to authenticate or fetch manually."; exit 1;
    }
  fi
  unzip -o GeminiMol.zip -d ./
  # Normalize layout: if a nested 'GeminiMol/' folder exists, move its contents up
  if [ -d GeminiMol ]; then
    echo "Normalizing directory layout: flattening nested GeminiMol/"
    shopt -s dotglob || true
    mv -f GeminiMol/* ./ || true
    shopt -u dotglob || true
    rmdir GeminiMol || true
  fi
  echo "✅ Unzipped GeminiMol.zip and normalized layout in $(pwd)"
elif [ $have_hf -eq 1 ]; then
  echo "Using huggingface-cli to fetch from AlphaMWang/GeminiMol"
  echo "Note: you may need to run 'huggingface-cli login' first."
  huggingface-cli download AlphaMWang/GeminiMol --repo-type model --local-dir ./ --local-dir-use-symlinks False || {
    echo "huggingface-cli download failed. Try 'huggingface-cli login' or use gdown."; exit 1;
  }
  echo "✅ Downloaded model files via Hugging Face"
else
  cat <<EOM
No downloader detected.
Please either:
  - Install gdown:  pip install gdown
    Then run this script again; or
  - Install Hugging Face CLI: pip install huggingface_hub
    huggingface-cli login
    Then run this script again; or
  - Manually download GeminiMol.zip (see README) to $(pwd) and unzip it here.
EOM
  exit 2
fi

echo "Done. Contents:"
ls -la
