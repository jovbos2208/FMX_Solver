#!/usr/bin/env bash
set -euo pipefail

# Fetch NRLMSIS2 and HWM14 sources into atm/extern/
# Usage:
#   NRLMSIS_URL=<tar.gz or git clone URL> \
#   HWM14_URL=<tar.gz or git clone URL> \
#   ./scripts/fetch_models.sh
#
# Notes:
# - This script only downloads/unpacks into atm/extern/. It does not build.
# - Choose official sources for licensing correctness.
# - For git URLs, ensure git is installed; for tar.gz, ensure curl/tar are available.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EXTERN_DIR="$ROOT_DIR/atm/extern"
mkdir -p "$EXTERN_DIR"

fetch_git() {
  local url="$1" dest="$2"
  if [[ -d "$dest/.git" ]]; then
    echo "[fetch] Updating $dest from $url";
    git -C "$dest" fetch --all && git -C "$dest" pull --ff-only
  else
    rm -rf "$dest"
    echo "[fetch] Cloning $url into $dest";
    git clone --depth=1 "$url" "$dest"
  fi
}

fetch_targz() {
  local url="$1" dest="$2"
  local tmp
  tmp="$(mktemp -d)"
  echo "[fetch] Downloading $url"
  curl -fsSL "$url" -o "$tmp/src.tgz"
  mkdir -p "$dest" && rm -rf "$dest"/*
  tar -xzf "$tmp/src.tgz" -C "$dest" --strip-components=1 || tar -xzf "$tmp/src.tgz" -C "$dest"
  rm -rf "$tmp"
}

if [[ -n "${NRLMSIS_URL:-}" ]]; then
  DEST="$EXTERN_DIR/nrlmsis2"
  if [[ "$NRLMSIS_URL" == *.git ]]; then
    fetch_git "$NRLMSIS_URL" "$DEST"
  else
    fetch_targz "$NRLMSIS_URL" "$DEST"
  fi
  echo "[ok] NRLMSIS2 fetched to $DEST"
else
  echo "[skip] Set NRLMSIS_URL to fetch NRLMSIS2 (git or tar.gz URL)."
fi

if [[ -n "${HWM14_URL:-}" ]]; then
  DEST="$EXTERN_DIR/hwm14"
  if [[ "$HWM14_URL" == *.git ]]; then
    fetch_git "$HWM14_URL" "$DEST"
  else
    fetch_targz "$HWM14_URL" "$DEST"
  fi
  echo "[ok] HWM14 fetched to $DEST"
else
  echo "[skip] Set HWM14_URL to fetch HWM14 (git or tar.gz URL)."
fi

echo "Done. You can now wire these in CMake or wrappers."

