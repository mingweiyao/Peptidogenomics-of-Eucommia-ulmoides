#!/bin/bash

# 增强版RNA-seq自动分析流程
# 作者：姚明伟
# 版本：2.0
# 更新：集成精准读长检测和错误处理机制

WORK_BASE_DIR="/data/Eu/NCBI/bam"
RAW_DATA_DIR="/media/wanglab/wanglab/EuRNASeq_NCBI/Transcriptomic/raw"
GENOME_DIR="/data/Eu/NCBI/Eu_genome"
GENOME_FILE="GWHBISF00000000.genome.fasta"
GTF_FILE="Eu.gtf"
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
    local index_dir="$GENOME_DIR/star_index_${sjdb_overhang}"
    
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
        STAR --runThreadN $THREADS \
             --genomeDir "$ref_index" \
             --readFilesIn "$input1" "$input2" \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix "/data/Eu/NCBI/bam/02alignment/${base}_" \
             --quantMode TranscriptomeSAM \
             --outSAMattributes Standard \
             --outSAMunmapped Within \
             --sjdbOverhang "$sjdb_overhang"
    else
        echo "$(timestamp) Aligning single-end file: $input1"
        STAR --runThreadN $THREADS \
             --genomeDir "$ref_index" \
             --readFilesIn "$input1" \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix "/data/Eu/NCBI/bam/02alignment/${base}_" \
             --quantMode TranscriptomeSAM \
             --outSAMattributes Standard \
             --outSAMunmapped Within \
             --sjdbOverhang "$sjdb_overhang" || return 1
    fi

    samtools index "/data/Eu/NCBI/bam/02alignment/${base}_Aligned.sortedByCoord.out.bam"
    local current_bam="/data/Eu/NCBI/bam/02alignment/${base}_Aligned.sortedByCoord.out.bam"
    local current_transcriptome_bam="/data/Eu/NCBI/bam/02alignment/${base}_Aligned.toTranscriptome.out.bam"
    if [[ -n "$input2" ]]; then
    	rsem-calculate-expression --bam --no-bam-output --paired-end -p $THREADS -q "$current_transcriptome_bam" "$GENOME_DIR/rsem_index/rsem_index" "/data/Eu/NCBI/bam/rsem/$base"
    	featureCounts -T 64 -p -t exon -g gene_id -a "$GENOME_DIR/$GTF_FILE" -o "/data/Eu/NCBI/bam/fcount/${base}_counts.txt" "$current_bam"
    else
    	rsem-calculate-expression --bam --no-bam-output -p $THREADS -q "$current_transcriptome_bam" "$GENOME_DIR/rsem_index/rsem_index" "/data/Eu/NCBI/bam/rsem/$base"
    	featureCounts -T 64 -t exon -g gene_id -a "$GENOME_DIR/$GTF_FILE" -o "/data/Eu/NCBI/bam/fcount/${base}_counts.txt" "$current_bam"
    fi
}

process_single_project() {
    local FILE_NAME="$1"
    local base=$(basename "$FILE_NAME" .sra)
    
    # 步骤1：解压原始数据
    echo "$(timestamp) Unpacking raw_data..."
    fastq-dump --split-3 -O "${WORK_BASE_DIR}/00raw" "$RAW_DATA_DIR/$FILE_NAME" || return 1
    
    # 步骤2：质控
    echo "$(timestamp) Trimming and filtering reads..."
    local is_paired=0
    if [[ -f "${WORK_BASE_DIR}/00raw/${base}_2.fastq" ]]; then
        is_paired=1
        fastp -i "${WORK_BASE_DIR}/00raw/${base}_1.fastq" \
              -I "${WORK_BASE_DIR}/00raw/${base}_2.fastq" \
              -o "${WORK_BASE_DIR}/01clean/${base}_1_trimmed.fastq" \
              -O "${WORK_BASE_DIR}/01clean/${base}_2_trimmed.fastq" \
              -w $((THREADS/2)) \
              --n_base_limit 5 \
              --qualified_quality_phred 15 \
              --unqualified_percent_limit 50 \
              --length_required 30 \
              --detect_adapter_for_pe \
              --dedup \
              --json "${WORK_BASE_DIR}/01clean/${base}_fastp.json" \
              --html "${WORK_BASE_DIR}/01clean/${base}_fastp.html" || return 1
    else
        fastp -i "${WORK_BASE_DIR}/00raw/${base}.fastq" \
              -o "${WORK_BASE_DIR}/01clean/${base}_trimmed.fastq" \
              -w $((THREADS/2)) \
              --n_base_limit 5 \
              --qualified_quality_phred 15 \
              --unqualified_percent_limit 50 \
              --length_required 30 \
              --dedup \
              --json "${WORK_BASE_DIR}/01clean/${base}_fastp.json" \
              --html "${WORK_BASE_DIR}/01clean/${base}_fastp.html" || return 1
    fi
    
    # 步骤3：读长检测和比对
    echo "$(timestamp) Determining read length..."
    local is_paired=0
    if [[ -f "${WORK_BASE_DIR}/01clean/${base}_1_trimmed.fastq" ]]; then
    	is_paired=1
    	input1="${WORK_BASE_DIR}/01clean/${base}_1_trimmed.fastq"
    	input2="${WORK_BASE_DIR}/01clean/${base}_2_trimmed.fastq"
    else
    	input1="${WORK_BASE_DIR}/01clean/${base}_trimmed.fastq"
    	input2=""
    fi
    local sjdb_overhang=$(get_sjdb_overhang "$input1")
    echo "$(timestamp) Using sjdbOverhang: $sjdb_overhang"
    
    # 步骤4：STAR比对
    local ref_index=$(get_or_build_star_index "$sjdb_overhang") || return 1
    run_star_alignment "$input1" "$input2" "$base" "$ref_index" "$sjdb_overhang" || return 1
    
    echo "$(timestamp) Successfully processed sample: $base"
    return 0
}

main() {
	for file in $(ls "$RAW_DATA_DIR"); do
        echo "Processing project: $file"
        process_single_project "$file"
    done
}

main
