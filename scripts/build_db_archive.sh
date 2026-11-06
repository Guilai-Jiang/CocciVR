#!/usr/bin/env bash
# Bundle the CocciVR reference database into a release asset.
#
#   ./scripts/build_db_archive.sh [DB_DIR] [OUTPUT.tar.gz]
set -euo pipefail

DB_DIR="${1:-db}"
OUT="${2:-coccivr-db-1.0.tar.gz}"

FILES=(
  VFDB_setB_pro.fas
  whole_cassette_SCCmec_database_EXTENDED_20171117.fasta
)

for f in "${FILES[@]}"; do
  if [[ ! -f "$DB_DIR/$f" ]]; then
    echo "Missing $DB_DIR/$f" >&2
    exit 1
  fi
done

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

for f in "${FILES[@]}"; do
  cp "$DB_DIR/$f" "$TMP/"
done

# Prebuilt BLAST indices for the SCCmec database, if present.
for ext in nhr nin nsq; do
  src="$DB_DIR/whole_cassette_SCCmec_database_EXTENDED_20171117.fasta.$ext"
  [[ -f "$src" ]] && cp "$src" "$TMP/"
done

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cp "$ROOT/docs/DATABASE-LICENSES.md" "$TMP/"

tar -czf "$OUT" -C "$TMP" .
echo "Wrote $OUT"
sha256sum "$OUT"
