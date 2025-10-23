#!/bin/bash

WORK_BASE_DIR="/data/Eu/Eu_rnaseq"
RAW_DATA_DIR="/data/Eu/Eu_rnaseq/sra"
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
    local index_dir="$GENOME_DIR/star_index_${sjdb_overhang}"
    
    if [[ -d "$index_dir" && -f "$index_dir/Genome" ]]; then
        echo "$index_dir"
        return 0
    fi
    
    echo "$(timestamp) Building new STAR index with sjdbOverhang=${sjdb_overhang}..."
    mkdir -p "$index_dir" || return 1
    
    /home/wanglab/anaconda3/envs/rnaseq/bin/STAR --runThreadN $THREADS \
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
        /home/wanglab/anaconda3/envs/rnaseq/bin/STAR --runThreadN $THREADS \
             --genomeDir "$ref_index" \
             --readFilesIn "$input1" "$input2" \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix "${WORK_BASE_DIR}/03alignment/${base}_" \
             --quantMode TranscriptomeSAM \
             --outSAMattributes Standard \
             --outSAMunmapped Within \
             --sjdbOverhang "$sjdb_overhang"
    else
        echo "$(timestamp) Aligning single-end file: $input1"
        /home/wanglab/anaconda3/envs/rnaseq/bin/STAR --runThreadN $THREADS \
             --genomeDir "$ref_index" \
             --readFilesIn "$input1" \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix "${WORK_BASE_DIR}/03alignment/${base}_" \
             --quantMode TranscriptomeSAM \
             --outSAMattributes Standard \
             --outSAMunmapped Within \
             --sjdbOverhang "$sjdb_overhang" || return 1
    fi

    samtools index "${WORK_BASE_DIR}/03alignment/${base}_Aligned.sortedByCoord.out.bam"
}

process_single_project() {
    local FILE_NAME="$1"
    local base=$(basename "$FILE_NAME" .sra)
    
    # # 步骤2: 质控
    # echo "$(timestamp) Generating quality reports for raw data..."
    # if [[ -f "${WORK_BASE_DIR}/00raw/${base}_2.fastq" ]]; then
    #     fastqc -t "$THREADS" -o "${WORK_BASE_DIR}/01quality/reports_before" \
    #             "${WORK_BASE_DIR}/00raw/${base}_1.fastq" \
    #             "${WORK_BASE_DIR}/00raw/${base}_2.fastq"
    # else
    #     fastqc -t "$THREADS" -o "${WORK_BASE_DIR}/01quality/reports_before" \
    #             "${WORK_BASE_DIR}/00raw/${base}.fastq" 
    # fi

    # # 步骤3：去接头
    # echo "$(timestamp) Trimming adapters..."
    # if [[ -f "${WORK_BASE_DIR}/00raw/${base}_2.fastq" ]]; then
    #     /data/Eu/fastp -i "${WORK_BASE_DIR}/00raw/${base}_1.fastq" \
    #           -I "${WORK_BASE_DIR}/00raw/${base}_2.fastq" \
    #           -o "${WORK_BASE_DIR}/02clean/${base}_1_trimmed.fastq" \
    #           -O "${WORK_BASE_DIR}/02clean/${base}_2_trimmed.fastq" \
    #           -w $((THREADS/2)) \
    #           --n_base_limit 5 \
    #           --qualified_quality_phred 15 \
    #           --unqualified_percent_limit 50 \
    #           --length_required 30 \
    #           --detect_adapter_for_pe \
    #           --dedup \
    #           --json "${WORK_BASE_DIR}/02clean/read_file/${base}_fastp.json" \
    #           --html "${WORK_BASE_DIR}/02clean/read_file/${base}_fastp.html" || return 1
    # else
    #     /data/Eu/fastp -i "${WORK_BASE_DIR}/00raw/${base}.fastq" \
    #           -o "${WORK_BASE_DIR}/02clean/${base}_trimmed.fastq" \
    #           -w $((THREADS/2)) \
    #           --n_base_limit 5 \
    #           --qualified_quality_phred 15 \
    #           --unqualified_percent_limit 50 \
    #           --length_required 30 \
    #           --dedup \
    #           --json "${WORK_BASE_DIR}/02clean/read_file/${base}_fastp.json" \
    #           --html "${WORK_BASE_DIR}/02clean/read_file/${base}_fastp.html" || return 1
    # fi
    
    # # 步骤4：去接头后质控
    # echo "$(timestamp) Generating quality reports for trimmed data..."
    # if [[ -f "${WORK_BASE_DIR}/02clean/${base}_2_trimmed.fastq" ]]; then
    #     fastqc -t "$THREADS" -o "${WORK_BASE_DIR}/01quality/reports_after" \
    #             "${WORK_BASE_DIR}/02clean/${base}_1_trimmed.fastq" \
    #             "${WORK_BASE_DIR}/02clean/${base}_2_trimmed.fastq"
    # else
    #     fastqc -t "$THREADS" -o "${WORK_BASE_DIR}/01quality/reports_after" \
    #             "${WORK_BASE_DIR}/02clean/${base}_trimmed.fastq"
    # fi

    # # 步骤5：读长检测和比对
    # echo "$(timestamp) Determining read length..."
    # if [[ -f "${WORK_BASE_DIR}/02clean/${base}_1_trimmed.fastq" ]]; then
    #     input1="${WORK_BASE_DIR}/02clean/${base}_1_trimmed.fastq"
    #     input2="${WORK_BASE_DIR}/02clean/${base}_2_trimmed.fastq"
    # else
    #     input1="${WORK_BASE_DIR}/02clean/${base}_trimmed.fastq"
    #     input2=""
    # fi

    # local sjdb_overhang=$(get_sjdb_overhang "$input1")
    # echo "$(timestamp) Using sjdbOverhang: $sjdb_overhang"
    
    # # 步骤6：STAR比对
    # local ref_index=$(get_or_build_star_index "$sjdb_overhang") || return 1
    # run_star_alignment "$input1" "$input2" "$base" "$ref_index" "$sjdb_overhang" || return 1
    
    # # 步骤7: 注释基因定量
    # echo "$(timestamp) Quantitative of annotated genes"
    # local current_bam="${WORK_BASE_DIR}/03alignment/${base}_Aligned.sortedByCoord.out.bam"
    # local current_transcriptome_bam="${WORK_BASE_DIR}/03alignment/${base}_Aligned.toTranscriptome.out.bam"
    # if [[ -n "$input2" ]]; then
    #     rsem-calculate-expression --bam --no-bam-output --paired-end -p $THREADS -q "$current_transcriptome_bam" "$GENOME_DIR/rsem_index/rsem_index" "${WORK_BASE_DIR}/04quan_annotated_gene/rsem/$base"
    #     featureCounts -T 64 -p -t exon -g gene_id -a "$GENOME_DIR/$GTF_FILE" -o "${WORK_BASE_DIR}/04quan_annotated_gene/fcount/${base}_counts.txt" "$current_bam"
    # else
    #     rsem-calculate-expression --bam --no-bam-output -p $THREADS -q "$current_transcriptome_bam" "$GENOME_DIR/rsem_index/rsem_index" "${WORK_BASE_DIR}/04quan_annotated_gene/rsem/$base"
    #     featureCounts -T 64 -t exon -g gene_id -a "$GENOME_DIR/$GTF_FILE" -o "${WORK_BASE_DIR}/04quan_annotated_gene/fcount/${base}_counts.txt" "$current_bam"
    # fi

    # # 步骤8: 新转录本鉴定
    # echo "$(timestamp) New transcript identification..."
    # local current_bam="${WORK_BASE_DIR}/03alignment/${base}_Aligned.sortedByCoord.out.bam"

    # stringtie -p "$THREADS" \
    # -o "${WORK_BASE_DIR}/05new_transcript/${base}.gtf" \
    # -G "$GENOME_DIR/$GTF_FILE" -B -l "STRG_${base}" \
    # "$current_bam"

    # # 提取转录本序列
    local output_dir="${WORK_BASE_DIR}/06new_transcript_positive"

    # gffread "${WORK_BASE_DIR}/05new_transcript/${base}.gtf" \
    # -g "${GENOME_DIR}/${GENOME_FILE}" \
    # -w "${output_dir}/${base}_transcripts.fasta" -W

    # echo "$(timestamp) Predicting coding potential for novel transcripts..."
    # # 1. CPC2 - 编码潜能评估
    # echo "$(timestamp) Running CPC2..."
    # python /home/wanglab/software/CPC2/bin/CPC2.py \
    # -i "${output_dir}/${base}_transcripts.fasta" \
    # -o "${output_dir}/cpc2_results/${base}_cpc2.txt"
    # # 2. PLEK - k-mer 方法
    # echo "$(timestamp) Running PLEK..."
    # python /home/wanglab/software/PLEK/PLEK.py \
    # -fasta "${output_dir}/${base}_transcripts.fasta" \
    # -out "${output_dir}/plek_results/${base}_plek.txt" \
    # -thread $THREADS

    # echo "$(timestamp) Running TransDecoder for sample: $base"
    cd "$output_dir"
    # TransDecoder.LongOrfs -t "${base}_transcripts.fasta" -m 50 -O "${base}_td"
    # TransDecoder.Predict  -t "${base}_transcripts.fasta" -O "${base}_td" --single_best_orf

    gtf_to_alignment_gff3.pl "${WORK_BASE_DIR}/05new_transcript/${base}.gtf" \
    > "${output_dir}/${base}.gff3"

    cdna_alignment_orf_to_genome_orf.pl \
    "${output_dir}/${base}_td/${base}_transcripts.fasta.transdecoder.gff3" \
    "${output_dir}/${base}.gff3" \
    "${output_dir}/${base}_transcripts.fasta" \
    > "${output_dir}/${base}.transdecoder.genome.gff3"

    gffread "${output_dir}/${base}.transdecoder.genome.gff3" \
    -g "${GENOME_DIR}/${GENOME_FILE}" \
    -x "${output_dir}/${base}.transdecoder.genome.cds.fa" \
    -y "${output_dir}/${base}.transdecoder.genome.pep.fa"

    echo "$(timestamp) Successfully processed sample: $base"
    return 0
}

main() {
    # 步骤0: 创建文件夹
    # echo "$(timestamp) Building dir..."
    # mkdir -p "${WORK_BASE_DIR}"/{00raw,01quality/{reports_before,reports_after},02clean/read_file,03alignment,04quan_annotated_gene/{rsem,fcount},05new_transcript,06new_transcript_positive/{cpc2_results,cnci_results,plek_results}}
    
    for file in $(ls "$RAW_DATA_DIR"); do
        # 步骤1：解压原始数据
    	echo "$(timestamp) Unpacking raw_data..."
    	fastq-dump --split-3 -O "${WORK_BASE_DIR}/00raw" "$RAW_DATA_DIR/$file" || return 1
    done
    
    echo "=== DEBUG INFO ==="
    echo "Script name: $0"
    echo "All arguments: $@"
    echo "Number of arguments: $#"
    echo "=== STARTING SCRIPT ==="
    # 清理可能冲突的变量
    unset STAR
    unset STAR-avx2  
    unset STAR-sse4.1
    
    for file in $(ls "$RAW_DATA_DIR"); do            
        echo "Processing project: $file"
        process_single_project "$file"
    done
}

main
