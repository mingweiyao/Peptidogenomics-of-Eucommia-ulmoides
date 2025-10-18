import os
import pandas as pd
from typing import List, Dict, Set
import argparse

def save_statistics(unique_df, gtf_files, combined_df, output_file):
    """保存统计信息"""
    stats_file = output_file.replace('.tsv', '_stats.txt')
    
    with open(stats_file, 'w') as f:
        f.write("转录本统计信息\n")
        f.write("=" * 50 + "\n")
        f.write(f"输入GTF文件数: {len(gtf_files)}\n")
        f.write(f"合并后总转录本数: {len(combined_df)}\n")
        f.write(f"去重后唯一转录本数: {len(unique_df)}\n")
        f.write(f"去重率: {(1 - len(unique_df)/len(combined_df))*100:.2f}%\n\n")
        
        # 源文件分布统计
        f.write("源文件分布统计:\n")
        file_counts = unique_df['source_file_count'].value_counts().sort_index()
        for count, num_transcripts in file_counts.items():
            f.write(f"  在 {count} 个文件中出现的转录本: {num_transcripts} 个\n")
        
        # 按染色体统计
        f.write("\n按染色体分布:\n")
        chrom_stats = unique_df['chromosome'].value_counts()
        for chrom, count in chrom_stats.items():
            f.write(f"  {chrom}: {count} 个转录本\n")
        
        # 按链统计
        f.write("\n按链分布:\n")
        strand_stats = unique_df['strand'].value_counts()
        for strand, count in strand_stats.items():
            f.write(f"  {strand}: {count} 个转录本\n")
        
        # 每个源文件中的转录本数量
        f.write("\n各源文件中的唯一转录本数量:\n")
        all_source_files = set()
        for files_str in unique_df['source_files']:
            all_source_files.update(files_str.split(', '))
        
        for source_file in sorted(all_source_files):
            count = unique_df['source_files'].str.contains(source_file).sum()
            f.write(f"  {source_file}: {count} 个唯一转录本\n")

def parse_gtf_file(file_path):
    """
    解析GTF文件，提取转录本信息
    """
    transcripts = []
    
    with open(file_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            # 跳过注释行
            if line.startswith('#'):
                continue
                
            parts = line.strip().split('\t')
            if len(parts) < 9:
                print(f"警告: 文件 {file_path} 第 {line_num} 行格式不正确，已跳过")
                continue
                
            chromosome, source, feature, start, end, score, strand, frame, attributes = parts
            
            # 只处理转录本相关的行（transcript, mRNA等）
            if feature in ['transcript', 'exon']:
                transcripts.append({
                    'chromosome': chromosome,
                    'start': int(start),
                    'end': int(end),
                    'strand': strand,
                    'source_file': os.path.basename(file_path)
                })
    
    return pd.DataFrame(transcripts)

def get_unique_transcripts_with_sources(gtf_folder, output_file):
    gtf_files = []
    for file in os.listdir(gtf_folder):
        if file.endswith('.gtf'):
            gtf_files.append(os.path.join(gtf_folder, file))
    
    print(f"找到 {len(gtf_files)} 个GTF文件")
    
    # 读取所有GTF文件并收集源文件信息
    all_transcripts = []
    source_records = {}  # 用于记录每个位置对应的源文件
    
    for gtf_file in gtf_files:
        try:
            df = parse_gtf_file(gtf_file)
            if not df.empty:
                all_transcripts.append(df)
                
                # 记录每个转录本位置对应的源文件
                for _, row in df.iterrows():
                    key = (row['chromosome'], row['start'], row['end'], row['strand'])
                    if key not in source_records:
                        source_records[key] = set()
                    source_records[key].add(os.path.basename(gtf_file))
                    
        except Exception as e:
            print(f"处理文件 {gtf_file} 时出错: {e}")
    
    if not all_transcripts:
        print("未成功提取到任何转录本信息")
        return pd.DataFrame()
    
    # 合并所有数据
    combined_df = pd.concat(all_transcripts, ignore_index=True)
    print(f"合并后总转录本数: {len(combined_df)}")
    
    # 根据染色体、起始位置、终止位置和链信息去重
    unique_df = combined_df.drop_duplicates(
        subset=['chromosome', 'start', 'end', 'strand'], 
        keep='first'
    ).copy()
    
    # 添加源文件列表列
    unique_df['source_files'] = unique_df.apply(
        lambda row: ', '.join(sorted(source_records.get(
            (row['chromosome'], row['start'], row['end'], row['strand']), set()
        ))), 
        axis=1
    )
    
    # 添加源文件数量列
    unique_df['source_file_count'] = unique_df.apply(
        lambda row: len(source_records.get(
            (row['chromosome'], row['start'], row['end'], row['strand']), set()
        )), 
        axis=1
    )
    
    # 排序
    unique_df = unique_df.sort_values(['chromosome', 'start', 'end', "strand"]).reset_index(drop=True)
    
    print(f"去重后唯一转录本数: {len(unique_df)}")
    
    # 保存结果
    if output_file:
        # 确保输出目录存在
        os.makedirs(os.path.dirname(output_file) if os.path.dirname(output_file) else '.', exist_ok=True)
        
        # 保存为TSV文件
        columns_order = ['chromosome', 'start', 'end', 'strand', 'transcript_id', 
                        'source_file_count', 'source_files', 'source_file']
        unique_df.to_csv(output_file, sep='\t', index=False, columns=columns_order)
        print(f"结果已保存到: {output_file}")
        
        # 同时保存统计信息
        save_statistics(unique_df, gtf_files, combined_df, output_file)
    
    return unique_df

def main():
    input_foder = r"G:\test_fractionation\20251018 imeta\file\00raw\new_transcript"
    output_file = r"D:\Desktop\Code\Peptidogenomics-of-Eucommia-ulmoides\00raw\merged_unique_transcripts.tsv"
    get_unique_transcripts_with_sources(input_foder, output_file)

if __name__ == "__main__":
    main()