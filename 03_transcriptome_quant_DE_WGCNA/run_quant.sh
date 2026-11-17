#!/bin/bash
set -uo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
BASE="${EU_RNASEQ_DIR:-$REPO_ROOT/data/rnaseq}"
STRINGTIE="${STRINGTIE:-stringtie}"
GFFCOMPARE="${GFFCOMPARE:-gffcompare}"

MERGED_GTF=$BASE/05new_transcript/merged.gff3
REF_GTF="${REF_GTF:-$REPO_ROOT/00raw_data/GWHBISF00000000.gtf}"
GTF_LIST=$BASE/05new_transcript/gtf_list.txt
BAM_DIR=$BASE/03alignment_gene/data_trans

OUT=$BASE/07new_transcript_quant
QUANT_DIR=$OUT/quant
LOG_DIR=$OUT/logs
CLASS_DIR=$OUT/gffcompare

MAX_PARALLEL=${MAX_PARALLEL:-10}
PER_THREAD=${PER_THREAD:-10}
LIMIT=${LIMIT:-0}

ts() { date '+%Y-%m-%d %H:%M:%S'; }
mkdir -p "$QUANT_DIR" "$LOG_DIR" "$CLASS_DIR"


sed 's/_transcript\.gtf$//' "$GTF_LIST" | sort -u > "$OUT/samples.all.txt"
: > "$OUT/samples.valid.txt"
: > "$OUT/missing_bams.txt"
while IFS= read -r s; do
    [ -z "$s" ] && continue
    if [ -s "$BAM_DIR/${s}_Aligned.sortedByCoord.out.bam" ]; then
        echo "$s" >> "$OUT/samples.valid.txt"
    else
        echo "$s" >> "$OUT/missing_bams.txt"
    fi
done < "$OUT/samples.all.txt"
echo "[$(ts)] Total samples: $(wc -l < "$OUT/samples.all.txt"); valid: $(wc -l < "$OUT/samples.valid.txt"); missing BAM: $(wc -l < "$OUT/missing_bams.txt")"


if [ ! -s "$CLASS_DIR/gffcmp.tracking" ]; then
    echo "[$(ts)] Stage1 running gffcompare..."
    "$GFFCOMPARE" -r "$REF_GTF" -o "$CLASS_DIR/gffcmp" "$MERGED_GTF" \
        > "$LOG_DIR/gffcompare.log" 2>&1
    echo "[$(ts)] Stage1 completed"
else
    echo "[$(ts)] Stage1 output already exists; skipping"
fi


run_one() {
    local s=$1
    local bam="$BAM_DIR/${s}_Aligned.sortedByCoord.out.bam"
    local odir="$QUANT_DIR/$s"
    if [ -s "$odir/$s.gtf" ]; then
        echo "[$(ts)] [skip] $soutput already exists"
        return 0
    fi
    mkdir -p "$odir"
    if "$STRINGTIE" -p "$PER_THREAD" -e -G "$MERGED_GTF" \
            -o "$odir/$s.gtf" -b "$odir" "$bam" \
            > "$LOG_DIR/$s.stringtie.log" 2>&1; then
        echo "[$(ts)] [done] $s"
    else
        echo "[$(ts)] [FAIL] $s (see $LOG_DIR/$s.stringtie.log)"
        return 1
    fi
}

n_done=0; n_skip=0; n_fail=0
total=$(wc -l < "$OUT/samples.valid.txt")
if [ "$LIMIT" -gt 0 ]; then
    head -n "$LIMIT" "$OUT/samples.valid.txt" > "$OUT/samples.batch.txt"
else
    cp "$OUT/samples.valid.txt" "$OUT/samples.batch.txt"
fi

echo "[$(ts)] Stage2 starting: batch of $(wc -l < "$OUT/samples.batch.txt") samples; parallelism $MAX_PARALLEL x $PER_THREAD threads per sample"
while IFS= read -r s; do
    [ -z "$s" ] && continue
    while [ "$(jobs -rp | wc -l)" -ge "$MAX_PARALLEL" ]; do
        wait -n || n_fail=$((n_fail+1))
    done
    if [ -s "$QUANT_DIR/$s/$s.gtf" ]; then n_skip=$((n_skip+1)); fi
    run_one "$s" &
done < "$OUT/samples.batch.txt"
wait
n_done=$(ls "$QUANT_DIR" | while IFS= read -r d; do [ -s "$QUANT_DIR/$d/$d.gtf" ] && echo x; done | wc -l)
echo "[$(ts)] Stage2 completed: $n_done / $total (skipped $n_skip); failed samples in $LOG_DIR"


echo "[$(ts)] Stage3 aggregating expression matrix..."
python3 "$OUT/aggregate_matrix.py" 2>&1 | tee "$LOG_DIR/aggregate.log"
echo "[$(ts)] All stages completed; output in$OUT"
