# ============================================================================
# [Fixed copy 2026-09-11] This file is a fixed copy of 00raw_data/aa_Eu_peptidegenomic_pipeline.py.
# Only fix: pd.concat(axis=1) in batch_analyze_peptides aligned rows by index label (see the
# [FIX 2026-09-11] comment inside the function); everything else is line-identical to the original.
# Re-running the pipeline with this version yields property columns that are correct by
# construction (Annotated_Peptides = 391,125 rows, removed_n = 3,321 rows; artifact rows gone),
# ============================================================================
from Bio import SeqIO
import pandas as pd
import logging
import re
from tqdm import tqdm
from intervaltree import IntervalTree
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import os
from collections.abc import Mapping
from openpyxl.utils.dataframe import dataframe_to_rows
from collections import defaultdict

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('peptide_analysis.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# 2. Load database information
def loading_peptide_db_information(peptide_db):
    peptide_db_dict = {}
    # >P1 GWHBISF00000001 1 |+| 1/110 1/330
    try:
        for rec in SeqIO.parse(peptide_db, "fasta"):
            db_id = rec.id
            parts = rec.description.split(" ")
            chrom = parts[1]
            strand = parts[3].split("|")[1]
            frame = parts[2]
            if strand == "-":
                frame = str(int(frame) + 3)
            start = parts[-1].split("/")[0]
            end = parts[-1].split("/")[1]
            peptide_db_dict[db_id] = (chrom, strand, frame, start, end, rec.seq)
    except Exception as e:
        logger.error(f"Failed to load database: {str(e)}")
    return  peptide_db_dict
def loading_protein_db_information(protein_db):
    protein_db_dict = {}
    '''
    >GWHPBISF000003	mRNA=GWHTBISF000003	Gene=GWHGBISF000003	Position=GWHBISF00000485: 221065-221100, 231722-231962, 232040-232285, 
                                    233270-233326, 233451-233711, 234204-234256, 234891-235005, 235094-235348, 238238-239319: +	
                                    Frame=0	OriID=evm.model.Chr1.10	OriTrascriptID=evm.model.Chr1.10	transl_table=1	OriGeneID=evm.TU.Chr1.10	OriSeqID=Chr1
    '''
    chrom_aa_dict = defaultdict(lambda: defaultdict(list))
    try:
        for rec in SeqIO.parse(protein_db, "fasta"):
            positions = []
            position_list = []
            mapped_positions = []
            db_id = rec.id
            parts = rec.description.split("\t")
            chrom_part = parts[3].split("=")[1].split(":")
            chrom = chrom_part[0]
            strand = chrom_part[-1].split(' ')[-1]
            position = chrom_part[1].split(", ")
            for r in position:
                start, end = map(int, r.split('-'))
                positions.append((start, end))
            for start, end in positions:
                position_list.extend(range(start, end + 1))
            if strand == "+":
                for i, char in enumerate(rec.seq):
                    mapped_positions.append((char, position_list[i * 3:i * 3 + 3]))
                    codon_pos = position_list[i * 3:i * 3 + 3]
                    chrom_aa_dict[chrom][char].append({"pos": codon_pos, "protein_id": db_id})
            else:
                position_list = position_list[:-3][::-1]
                for i, char in enumerate(rec.seq):
                    mapped_positions.append((char, position_list[i * 3:i * 3 + 3]))
                    codon_pos = position_list[i * 3:i * 3 + 3]
                    chrom_aa_dict[chrom][char].append({"pos": codon_pos, "protein_id": db_id})
            protein_db_dict[db_id] = (chrom, strand, mapped_positions, rec.seq)
    except Exception as e:
        logger.error(f"Failed to load database: {str(e)}")
    chrom_aa_dict_final = {
        chrom: dict(aa_pos) 
        for chrom, aa_pos in chrom_aa_dict.items()
    }
    return protein_db_dict, chrom_aa_dict_final

# 3. Chromosome length statistics
def genome_chrom_length(genome_file):
    # >GWHBISF00000001	OriSeqID=F000542F	Len=330
    chrom_length_dict = {}
    try:
        for rec in SeqIO.parse(genome_file, "fasta"):
            match = re.search(r'Len=(\d+)', rec.description)
            if match:
                chrom_length_dict[rec.id] = int(match.group(1))
        logger.info(f"Obtained lengths for {len(chrom_length_dict)} chromosomes")
    except Exception as e:
        logger.error(f"Failed to get chromosome lengths: {str(e)}")
        raise
    return chrom_length_dict

# 4. Compute peptide genomic coordinates
def calculate_peptide_coordinates(peptides, peptide_db_info, chrom_length_dict):
    peptides = peptides.copy()
    peptides["start"] = None
    peptides["end"] = None
    peptides["strand"] = None
    peptides["chrom"] = None
    peptides["frame"] = None
    for index, row in tqdm(peptides.iterrows(), total=len(peptides), desc="Computing peptide coordinates"):
        try:
            id = row['accessions']
            seq = row["sequence"]
            if id not in peptide_db_info:
                continue
            db_chrom, db_strand, db_frame, db_start, db_end, db_seq = peptide_db_info[id]
            chrom_length = chrom_length_dict[db_chrom]
            pos = db_seq.find(seq)
            if pos == -1:
                continue
            if db_strand == "+":
                start = int(db_start) + pos * 3
                end = start + len(seq) * 3 - 1
            else:
                start_strand = int(db_start) + pos * 3
                end_strand = start_strand + len(seq) * 3 - 1
                end = chrom_length - start_strand + 1
                start = chrom_length - end_strand + 1
            peptides.at[index, 'start'] = start
            peptides.at[index, 'end'] = end
            peptides.at[index, 'strand'] = db_strand
            peptides.at[index, 'chrom'] = db_chrom
            peptides.at[index, 'frame'] = db_frame
        except Exception as e:
            logger.warning(f"Row {index} failed: {str(e)}")
            continue
    return peptides
def calculate_protein_coordinates(proteins, protein_db_info):
    proteins = proteins.copy()
    proteins["start"] = None
    proteins["end"] = None
    proteins["strand"] = None
    proteins["chrom"] = None
    proteins["frame"] = 0
    for index, row in tqdm(proteins.iterrows(), total=len(proteins), desc="Computing protein coordinates"):
        try:
            protein_id = row['accessions']
            seq = row["sequence"]
            if protein_id not in protein_db_info:
                continue
            db_chrom, db_strand, db_mapped_positions, db_seq = protein_db_info[protein_id]
            seq_len = len(seq)
            pos = db_seq.find(seq)
            if pos == -1 or (pos + seq_len) > len(db_mapped_positions):
                continue
            first_pos = db_mapped_positions[pos][1]
            last_pos = db_mapped_positions[pos + seq_len - 1][1]
            if db_strand == "+":
                start = first_pos[0]
                end = last_pos[-1]
            else:
                start = last_pos[-1]
                end = first_pos[0]
            proteins.at[index, 'start'] = start
            proteins.at[index, 'end'] = end
            proteins.at[index, 'strand'] = db_strand
            proteins.at[index, 'chrom'] = db_chrom   
        except Exception as e:
            logger.warning(f"Row {index} failed: {str(e)}")
            continue       
    return proteins

# 5. Parse GFF file
def region_calculate(mrnas, exons):
    required_cols = ['chrom', 'type', 'strand', 'start', 'end']
    if not all(col in mrnas.columns for col in required_cols):
        raise ValueError("mRNAs DataFrame is missing required columns")
    if not all(col in exons.columns for col in required_cols):
        raise ValueError("Exons DataFrame is missing required columns")
    results = []
    mrnas_group = mrnas.groupby(['chrom', 'strand'])
    exons_group = exons.groupby(['chrom', 'strand'])
    for (chrom, strand), mrna_group in mrnas_group:
        tree = IntervalTree()    
        for _, row in mrna_group.iterrows():
            tree.addi(row["start"], row["end"] + 1)
        if (chrom, strand) in exons_group.groups:
            for _, exon_row in exons_group.get_group((chrom, strand)).iterrows():
                tree.chop(exon_row["start"], exon_row["end"] + 1)
        for interval in sorted(tree):
            results.append({
                "chrom": chrom,
                "type": "intron",
                "start": interval.begin,
                "end": interval.end - 1,
                "strand": strand
            })   
    return pd.DataFrame(results)
def intergenic_calculate(chrom_length_dict, genes, strand):
    if not isinstance(chrom_length_dict, dict):
        raise TypeError("chrom_length_dict must be a dict")
    required_cols = ['chrom', 'start', 'end', 'strand']
    if not all(col in genes.columns for col in required_cols):
        raise ValueError(f"genes DataFrame is missing required columns: {required_cols}")
    if strand not in ('+', '-'):
        raise ValueError("strand must be '+' or '-'")
    intergenic_regions = []
    for chrom, length in chrom_length_dict.items():
        try:
            chrom_length = int(length)
            if chrom_length <= 0:
                raise ValueError(f"Chromosome {chrom} length must be a positive integer")
        except (ValueError, TypeError):
            logger.warning(f"Ignoring invalid chromosome length: {chrom}={length}")
            continue
        tree = IntervalTree()
        tree.addi(1, chrom_length + 1)
        chrom_genes = genes[(genes['chrom'] == chrom)]
        if not chrom_genes.empty:
            for _, gene in chrom_genes.iterrows():
                try:
                    tree.chop(gene["start"], gene["end"] + 1)
                except ValueError as e:
                    logger.warning(f"Ignoring invalid gene interval {chrom}:{gene['start']}-{gene['end']}: {str(e)}")
                    continue
        tree.merge_overlaps(strict=False)
        for interval in sorted(tree):
            intergenic_regions.append({
                "chrom": chrom,
                "type": "intergenic",
                "start": interval.begin,
                "end": interval.end - 1,
                "strand": strand
            })
    return pd.DataFrame(intergenic_regions)
def parse_gff_file(gff_file, chrom_length_dict):
    try:
        gff_data = []
        with open(gff_file, "r") as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 9:
                    logger.warning(f"Ignoring incomplete line (too few columns): line {line_num}")
                    continue
                try:
                    record = {
                        'chrom': parts[0],
                        'type': parts[2],
                        'start': int(parts[3]),
                        'end': int(parts[4]),
                        'strand': parts[6]
                    }
                    if record['start'] > record['end']:
                        logger.warning(f"Ignoring invalid coordinates (start>end): line {line_num}")
                        continue
                    gff_data.append(record)
                except (ValueError, IndexError) as e:
                    logger.warning(f"Ignoring malformed line: line {line_num}: {str(e)}")
                    continue
        gff_df = pd.DataFrame(gff_data)
        if gff_df.empty:
            logger.error("GFF file contains no valid records")
            return {}
        def filter_features(df, feature_type, strand=None):
            mask = df['type'] == feature_type
            if strand is not None:
                mask &= df['strand'] == strand
            return df[mask].copy()
        mrnas_pos = filter_features(gff_df, 'mRNA', '+')
        mrnas_neg = filter_features(gff_df, 'mRNA', '-')
        exons_pos = filter_features(gff_df, 'exon', '+')
        exons_neg = filter_features(gff_df, 'exon', '-')
        cds_pos = filter_features(gff_df, 'CDS', '+')
        cds_neg = filter_features(gff_df, 'CDS', '-')
        utr3_pos = filter_features(gff_df, 'three_prime_UTR', '+')
        utr3_neg = filter_features(gff_df, 'three_prime_UTR', '-')
        utr5_pos = filter_features(gff_df, 'five_prime_UTR', '+')
        utr5_neg = filter_features(gff_df, 'five_prime_UTR', '-')
        intron_pos = region_calculate(mrnas_pos, exons_pos)
        intron_neg = region_calculate(mrnas_neg, exons_neg)
        intergenic_pos = intergenic_calculate(chrom_length_dict, mrnas_pos, "+")
        intergenic_neg = intergenic_calculate(chrom_length_dict, mrnas_neg, "-")
        features = pd.concat([
            cds_pos, cds_neg,
            utr3_pos, utr3_neg,
            utr5_pos, utr5_neg,
            intron_pos, intron_neg,
            intergenic_pos, intergenic_neg
        ], ignore_index=True)
        feature_dict = {}
        for _, row in features.iterrows():
            key = (row['chrom'], row['strand'])
            feature_dict.setdefault(key, []).append(
                (row['start'], row['end'], row['type'])
            )
        logger.info(
            f"GFF parsed: processed {len(gff_df)} records, "
            f"features built for {len(feature_dict)} chromosome/strand combinations"
        )
        return feature_dict
    except Exception as e:
        logger.error(f"Failed to parse GFF: {str(e)}")
        raise

# 6. Filter N-region overlaps
def detect_n_regions(genome_file):
    n_regions = []
    try:
        for rec in SeqIO.parse(genome_file, "fasta"):
            chrom = rec.description.split()[0]
            seq = str(rec.seq).upper()
            n_start = None
            for i, base in enumerate(seq, 1):
                if base == 'N':
                    if n_start is None:
                        n_start = i
                else:
                    if n_start is not None:
                        n_regions.append([chrom, n_start, i-1])
                        n_start = None
            if n_start is not None:
                n_regions.append([chrom, n_start, len(seq)])
        n_df = pd.DataFrame(n_regions, columns=['chrom', 'start', 'end'])
        logger.info(f"Found {len(n_df)} N-regions")
        return n_df
    except Exception as e:
        logger.error(f"N-region detection failed: {str(e)}")
        raise
def filter_n_overlaps(peptide_df, n_regions):
    required_cols = ['chrom', 'start', 'end']
    if not all(col in peptide_df.columns for col in required_cols):
        raise ValueError("Peptide DataFrame is missing required columns")
    if not all(col in n_regions.columns for col in required_cols):
        raise ValueError("N-region DataFrame is missing required columns")
    from intervaltree import IntervalTree
    chrom_trees = {}
    for _, row in n_regions.iterrows():
        chrom = row['chrom']
        if chrom not in chrom_trees:
            chrom_trees[chrom] = IntervalTree()
        chrom_trees[chrom].addi(row['start'], row['end']+1)
    peptide_df = peptide_df.copy()
    peptide_df['overlaps_n'] = False
    for idx, row in tqdm(peptide_df.iterrows(), total=len(peptide_df), desc="Filtering N-region overlaps"):
        chrom = row['chrom']
        start = row['start']
        end = row['end']
        if pd.isna(start) or pd.isna(end) or chrom not in chrom_trees:
            continue
        if chrom_trees[chrom].overlaps(start, end+1):
            peptide_df.at[idx, 'overlaps_n'] = True
    filtered = peptide_df[~peptide_df['overlaps_n']].copy()
    removed = peptide_df[peptide_df['overlaps_n']].copy()
    filtered.drop(columns=['overlaps_n'], inplace=True)
    removed.drop(columns=['overlaps_n'], inplace=True)
    logger.info(f"Filter result: kept {len(filtered)}, removed {len(removed)}")
    return filtered, removed

# 7. Annotate peptide genomic positions
def annotate_peptide_loc(gff_info, peptide_file_coord):
    required_cols = ['chrom', 'start', 'end', 'strand']
    if not all(col in peptide_file_coord.columns for col in required_cols):
        raise ValueError("Peptide coordinate DataFrame is missing required columns")
    type_list = []
    for _, row in tqdm(peptide_file_coord.iterrows(), total=len(peptide_file_coord), desc="Annotating genomic positions"):
        chrom = row['chrom']
        start = row['start']
        end = row['end']
        strand = row['strand']
        search_key = (chrom, strand)
        current_type = None
        junction_found = False
        if search_key in gff_info:
            features = gff_info[search_key]
            for feat_start, feat_end, feat_type in features:
                if feat_start <= start < end <= feat_end:
                    current_type = feat_type
                    break
                elif (start < feat_start < end) or (start < feat_end < end):
                    junction_found = True
            if current_type is None and junction_found:
                current_type = "junction"
        type_list.append(current_type)
    peptide_file_coord = peptide_file_coord.copy()
    peptide_file_coord['type'] = type_list
    return peptide_file_coord
def annotate_protein_loc(protein_file_coord, protein_db):
    position_dict = {}
    for record in SeqIO.parse(protein_db, "fasta"):
        pos_info = record.description.split("\t")[3].split(": ")
        position_info = pos_info[-2]
        strand = pos_info[-1]
        ranges = position_info.split(',')
        boundaries = []
        for r in ranges:
            start, end = r.strip().split('-')
            start = int(start)
            end = int(end)
            boundaries.append((start, end))
        position_dict[record.id] = {
            'strand': strand,
            'boundaries': sorted(boundaries)
        }
    type_list = []
    first_boundary = []
    second_boundary = []
    for _, row in tqdm(protein_file_coord.iterrows(), total=len(protein_file_coord), desc="Annotating genomic positions"):
        start = row['start']
        end = row['end']
        accession = row['accessions']
        chrom_data = position_dict[accession]
        boundaries = chrom_data['boundaries']
        if start != None:
            type_list.append("CPs")
            flag = 0
            for i in range(len(boundaries)-1):
                if (start <= boundaries[i][1]) and (end >= boundaries[i+1][0]):
                    first_boundary.append(boundaries[i][1])
                    second_boundary.append(boundaries[i+1][0])
                    flag = 1
                    break
            if flag==0:
                first_boundary.append(None)
                second_boundary.append(None)
        else:
            type_list.append("error")
            first_boundary.append(None)
            second_boundary.append(None)
    protein_file_coord['type'] = type_list
    protein_file_coord['first_boundary'] = first_boundary
    protein_file_coord['second_boundary'] = second_boundary
    return protein_file_coord

# 8. Classify CDS peptides
def cds_peptide_analysis(annotated_peptides, protein_db_info, chrom_aa_dict):
    required_cols = ['type', 'accessions', 'sequence', 'chrom', 'start', 'end', 'strand']
    if not all(col in annotated_peptides.columns for col in required_cols):
        raise ValueError(f"Input DataFrame is missing required columns: {required_cols}")
    cds_types = pd.Series(index=annotated_peptides.index, dtype=object)
    for idx, row in tqdm(annotated_peptides.iterrows(),
                        total=len(annotated_peptides),
                        desc="Classifying CDS peptides"):
        if row['type'] != 'CDS':
            cds_types[idx] = None
            continue
        peptide_id = row['accessions']
        seq = row['sequence']
        chrom = row['chrom']
        start = row['start']
        end = row['end']
        strand = row['strand']
        if peptide_id in protein_db_info:
            db_chrom, db_strand, mapped_position, db_seq = protein_db_info[peptide_id]
            if strand != db_strand:
                cds_types[idx] = "CDS(strand mismatch)"
                continue
            index = db_seq.find(seq)
            if index == -1:
                cds_types[idx] = "CDS(out of frame)"
                continue
            if strand == "+":
                seq_start = mapped_position[index][1][0]
                seq_end = mapped_position[index + len(seq) - 1][1][-1]
            else:
                seq_start = mapped_position[index + len(seq) - 1][1][-1]
                seq_end = mapped_position[index][1][0]    
            if seq_start == start and seq_end == end:
                cds_types[idx] = "CDS(verified by protein_db)"
            else:
                cds_types[idx] = "CDS(out of frame)"
        elif chrom in chrom_aa_dict:
            first_aa = seq[0]
            matched = False
            if first_aa in chrom_aa_dict[chrom]:
                for pos_info in chrom_aa_dict[chrom][first_aa]:
                    pos_group = pos_info["pos"]
                    ref_pos = pos_group[0] if strand == '+' else pos_group[-1]
                    if ref_pos == start:
                        matched_protein_id = pos_info["protein_id"] 
                        if matched_protein_id in protein_db_info:
                            db_chrom, db_strand, mapped_position, db_seq = protein_db_info[matched_protein_id]
                            if strand != db_strand:
                                cds_types[idx] = "CDS(strand mismatch)"
                                matched = True
                                break
                            index = db_seq.find(seq)
                            if index == -1:
                                cds_types[idx] = "CDS(out of frame)"
                                matched = True
                                break
                            if strand == "+":
                                seq_start = mapped_position[index][1][0]
                                seq_end = mapped_position[index + len(seq) - 1][1][-1]
                            else:
                                seq_start = mapped_position[index + len(seq) - 1][1][-1]
                                seq_end = mapped_position[index][1][0]    
                                
                            if seq_start == start and seq_end == end:
                                cds_types[idx] = "CDS(verified by protein_db)"
                            else:
                                cds_types[idx] = "CDS(out of frame)"
                            matched = True
                            break
                
                if not matched:
                    cds_types[idx] = "CDS(out of frame)"
            else:
                cds_types[idx] = "CDS(out of frame)"
        else:
            cds_types[idx] = "CDS(out of frame)"
    result_df = annotated_peptides.copy()
    result_df['cds_type'] = cds_types
    return result_df

# 9. Compute peptide properties
def analyze_peptide_properties(sequence):
    if not sequence or not isinstance(sequence, str):
        return {
            "length": 0,
            "molecular_weight": None,
            "isoelectric_point": None,
            "gravy": None,
            "aromaticity": None,
            "instability_index": None
        }
    try:
        analysis = ProteinAnalysis(sequence)
        properties = {
            "length": len(sequence),
            "molecular_weight": analysis.molecular_weight(),
            "isoelectric_point": analysis.isoelectric_point(),
            "gravy": analysis.gravy(),
        }
        try:
            properties.update({
                "instability_index": analysis.instability_index()
            })
        except Exception as e:
            logger.warning(f"Some property computations failed: {str(e)}")
        return properties  
    except Exception as e:
        logger.error(f"Peptide analysis failed: {str(e)}")
        return {
            "length": len(sequence),
            "molecular_weight": None,
            "isoelectric_point": None,
            "gravy": None,
            "aromaticity": None,
            "instability_index": None
        }
def batch_analyze_peptides(annotated_peptides_finally):
    if 'sequence' not in annotated_peptides_finally.columns:
        logger.error("Column 'sequence' not found in the Excel file")
        return None
    results = []
    for _, row in tqdm(annotated_peptides_finally.iterrows(), 
                          total=len(annotated_peptides_finally),
                          desc="Computing properties"):
        seq = row.get('sequence', '')
        results.append(analyze_peptide_properties(seq))
    result_df = pd.DataFrame(results)
    # [FIX 2026-09-11] Original code: pd.concat([annotated_peptides_finally, result_df], axis=1)
    # Bug: concat(axis=1) aligns on index labels. After filtering, cds_peptide_analysis /
    #   filter_n_overlaps keep the original labels (with gaps; removed even has duplicated labels),
    #   while result_df has a fresh RangeIndex, so label-based pairing is wrong and property
    #   columns attach to the wrong peptide rows (~92.3% of S2 rows misaligned). Duplicated
    #   labels also create artifact rows (Annotated_Peptides 394,431 physical vs 391,125 valid
    #   rows; removed_n 6,630 vs 3,321 - the gap is exactly these artifacts).
    base = annotated_peptides_finally.reset_index(drop=True)
    output_df = pd.concat([base, result_df], axis=1)
    return output_df

# 10. Save output files
def save_to_excel(data_dict, output_file):
    try:
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            for sheet_name, data in data_dict.items():
                if isinstance(data, Mapping):
                    if all(isinstance(v, (list, tuple)) for v in data.values()):
                        rows = []
                        for key, values in data.items():
                            chrom, strand = key
                            for (start, end, feat_type) in values:
                                rows.append([chrom, strand, start, end, feat_type])
                        df = pd.DataFrame(rows, columns=["chrom", "strand", "start", "end", "type"])
                    else:
                        df = pd.DataFrame.from_dict(data, orient='index')
                else:
                    df = data if isinstance(data, pd.DataFrame) else pd.DataFrame()
                if df.empty:
                    df = pd.DataFrame({"Info": [f"No data in {sheet_name}"]})
                safe_name = re.sub(r'[\\/*?\[\]:]', '_', str(sheet_name))[:31]
                df.to_excel(
                    writer,
                    sheet_name=safe_name,
                    index=False,
                    header=not df.empty,
                    freeze_panes=(1, 0)
                )
                if not df.empty:
                    worksheet = writer.sheets[safe_name]
                    for col in worksheet.columns:
                        max_len = max((
                            len(str(cell.value)) 
                            for cell in col
                        ), default=0) + 2
                        worksheet.column_dimensions[col[0].column_letter].width = min(max_len, 50)
        logger.info(f"Results saved to: {output_file}")
        return True
    except Exception as e:
        logger.error(f"Failed to save: {str(e)}")
        if os.path.exists(output_file):
            os.remove(output_file)
        return False

# Overall workflow
def run_workflow(peptide_db, protein_db, peptide_file, protein_file, gff_file, genome_file, output_dir):
    logger.info("=== Peptide analysis workflow started ===")
    try:
        logger.info("1. Loading peptide data...")
        peptides = pd.read_excel(peptide_file)
        peptides.columns = peptides.columns.str.strip().str.lower()
        proteins = pd.read_excel(protein_file)
        proteins.columns = proteins.columns.str.strip().str.lower()

        logger.info("2. Loading database data...")
        peptide_db_info = loading_peptide_db_information(peptide_db)
        protein_db_info, chrom_aa_dict = loading_protein_db_information(protein_db)

        logger.info("3. Chromosome length statistics")
        chrom_length_dict = genome_chrom_length(genome_file)

        logger.info("4. Computing peptide coordinates")
        peptide_file_coord = calculate_peptide_coordinates(peptides, peptide_db_info, chrom_length_dict)
        protein_file_coord = calculate_protein_coordinates(proteins, protein_db_info)

        logger.info("5. Parsing GFF file")
        gff_info = parse_gff_file(gff_file, chrom_length_dict)

        logger.info('6. Filtering N-region overlaps')
        n_regions = detect_n_regions(genome_file)
        filtered, removed = filter_n_overlaps(peptide_file_coord, n_regions)

        logger.info("7. Annotating peptide genomic positions")
        annotated_peptides = annotate_peptide_loc(gff_info, filtered)
        annotated_proteins = annotate_protein_loc(protein_file_coord, protein_db)

        logger.info("8. Classifying CDS peptides")
        annotated_peptides_finally = cds_peptide_analysis(annotated_peptides, protein_db_info, chrom_aa_dict)
        annotated_proteins_finally = cds_peptide_analysis(annotated_proteins, protein_db_info, chrom_aa_dict)

        logger.info("9. Computing peptide properties")
        peptide_finally = batch_analyze_peptides(annotated_peptides_finally)
        removed_n = batch_analyze_peptides(removed)
        protein_finally = batch_analyze_peptides(annotated_proteins_finally)

        logger.info("10. Saving results")
        output_data = {
            "Annotated_Peptides": peptide_finally,
            "Annotated_Proteins": protein_finally,
            "Genomic_Features": gff_info,
            "removed_n": removed_n
        }
        output_file = os.path.join(output_dir, "peptide_analysis_results.xlsx")
        save_to_excel(output_data, output_file)
        
        logger.info(f"=== Workflow completed ===")
    except Exception as e:
        logger.error(f"Workflow failed: {str(e)}")
        raise

if __name__ == "__main__":
    import os as _os
    RAW = _os.environ.get(
        "AA_RAW_DIR",
        r"./data/rejected_20260204/00raw_data")
    input_files = {
        "peptide_db": f"{RAW}/Eu_peptide_database_customized_5.fa",
        "protein_db": f"{RAW}/Eu_peptide_database.fa",
        "peptide_file": f"{RAW}/Eu_peptide_raw.xlsx",
        "protein_file": f"{RAW}/Eu_protein_raw.xlsx",
        "gff_file": f"{RAW}/GWHBISF00000000.gff",
        "genome_file": f"{RAW}/Eu_genome.fasta",
        "output_dir": _os.environ.get("AA_OUT_DIR", r"./output/_aa_pipeline_fixed")
    }
    try:
        run_workflow(**input_files)
    except Exception as e:
        logger.critical(f"Terminated: {str(e)}")