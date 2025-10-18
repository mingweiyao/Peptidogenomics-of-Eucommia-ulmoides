#!/bin/bash

WORK_BASE_DIR="/data/Eu/bam_for_sp/STAR_bam"
RAW_DATA_DIR="/data/Eu/test"
SRA_DATA_DIR="/media/wanglab/wanglab/Eu/EuRNASeq_NCBI/Transcriptomic/test"
GENOME_DIR="/data/Eu/Eu_genome"
GENOME_FILE="GWHBISF00000000.genome.fasta"
GTF_FILE="GWHBISF00000000.gtf"
THREADS=100

timestamp() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')]"
}

get_dominant_read_length() {
    local fastq=$1
    bioawk -c fastx '{
        len[length($seq)]++
    } END {
        max=0; mode=0;
        for (l in len) {
            if (len[l] > max) {
                max = len[l]
                mode = l
            }
        }
        print mode
    }' "$fastq"
}

get_sjdb_overhang() {
    local input_file="$1"
    local read_length=$(get_dominant_read_length "$input_file")
    local sjdb_overhang=$((read_length - 1))
    [[ $sjdb_overhang -lt 50 ]] && sjdb_overhang=50
    echo "$sjdb_overhang"
}

get_or_build_star_index() {
    local sjdb_overhang="$1"
    local index_dir="$GENOME_DIR/star_index_sp_${sjdb_overhang}"
    
    if [[ -d "$index_dir" && -f "$index_dir/Genome" ]]; then
        echo "$index_dir"
        return 0
    fi
    
    echo "$(timestamp) Building new STAR index with sjdbOverhang=${sjdb_overhang}..."
    mkdir -p "$index_dir" || return 1
    
    STAR --runThreadN $THREADS \
         --runMode genomeGenerate \
         --genomeDir "$index_dir" \
         --genomeFastaFiles "$GENOME_DIR/$GENOME_FILE" \
         --sjdbGTFfile "$GENOME_DIR/$GTF_FILE" \
         --sjdbOverhang $sjdb_overhang \
         --genomeSAindexNbases 13 || return 1
    
    echo "$index_dir"
}

run_star_alignment() {
    local input1=$1
    local input2=$2
    local base=$3
    local ref_index=$4
    local sjdb_overhang=$5

    if [[ -n "$input2" ]]; then
        echo "$(timestamp) Aligning paired-end files: $input1 and $input2"
        STAR-avx2 --runThreadN $THREADS \
             --genomeDir "$ref_index" \
             --readFilesIn "$input1" "$input2" \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix "/data/Eu/bam_for_sp/STAR_bam/${base}_" \
             --outSAMattributes Standard \
             --outFilterMultimapNmax 1 \
             --outSAMmultNmax 1 \
             --sjdbOverhang "$sjdb_overhang"
    else
        echo "$(timestamp) Aligning single-end file: $input1"
        STAR-avx2 --runThreadN $THREADS \
             --genomeDir "$ref_index" \
             --readFilesIn "$input1" \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix "/data/Eu/bam_for_sp/STAR_bam/${base}_" \
             --outSAMattributes Standard \
             --outFilterMultimapNmax 1 \
             --outSAMmultNmax 1 \
             --sjdbOverhang "$sjdb_overhang" || return 1
    fi

    samtools view -@ 50 -h -q 255 -b "/data/Eu/bam_for_sp/STAR_bam/${base}_Aligned.sortedByCoord.out.bam" > "/data/Eu/bam_for_sp/STAR_bam/${base}_unique_Aligned.sortedByCoord.out.bam"
    samtools index "/data/Eu/bam_for_sp/STAR_bam/${base}_unique_Aligned.sortedByCoord.out.bam"
    local current_bam="/data/Eu/bam_for_sp/STAR_bam/${base}_unique_Aligned.sortedByCoord.out.bam"
    if [[ -n "$input2" ]]; then
    	featureCounts -T 64 -Q 255 -p --countReadPairs -t exon -g gene_id -a "$GENOME_DIR/$GTF_FILE" -o "/data/Eu/bam_for_sp/STAR_bam/featurecount_gene/${base}_counts.txt" "$current_bam"
    else
    	featureCounts -T 64 -Q 255 -t exon -g gene_id -a "$GENOME_DIR/$GTF_FILE" -o "/data/Eu/bam_for_sp/STAR_bam/featurecount_gene/${base}_counts.txt" "$current_bam"
    fi
}

process_single_project() {
    local FILE_NAME="$1"
    local base=$(basename "$FILE_NAME" .sra)
    
    echo "$(timestamp) Determining read length..."
    if [[ -f "$RAW_DATA_DIR/${base}_1_trimmed.fastq" ]]; then
    	input1="$RAW_DATA_DIR/${base}_1_trimmed.fastq"
    	input2="$RAW_DATA_DIR/${base}_2_trimmed.fastq"
    else
    	input1="$RAW_DATA_DIR/${base}_trimmed.fastq"
    	input2=""
    fi
    local sjdb_overhang=$(get_sjdb_overhang "$input1")
    echo "$(timestamp) Using sjdbOverhang: $sjdb_overhang"
    
    local ref_index=$(get_or_build_star_index "$sjdb_overhang") || return 1
    run_star_alignment "$input1" "$input2" "$base" "$ref_index" "$sjdb_overhang" || return 1
    
    echo "$(timestamp) Successfully processed sample: $base"
    return 0
}

main() {
	for file in $(ls "$SRA_DATA_DIR"); do
        echo "Processing project: $file"
        process_single_project "$file"
    done
}

main

