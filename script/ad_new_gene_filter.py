import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from typing import Dict, List

def load_data(expression_matrix_path: str, sample_design_path: str) -> (pd.DataFrame, Dict[str, List[str]]):
    """
    加载表达矩阵和实验设计文件。
    """
    # 加载表达矩阵
    expression_df = pd.read_csv(expression_matrix_path, sep='\t', index_col=0)
    
    # 加载实验设计
    design_df = pd.read_csv(sample_design_path, sep='\t')
    
    # 构建重复组字典
    replicate_groups = {}
    for group_name in design_df['group'].unique():
        samples_in_group = design_df[design_df['group'] == group_name]['sample'].tolist()
        samples_in_group = [s for s in samples_in_group if s in expression_df.columns]
        if len(samples_in_group) > 1:
            replicate_groups[group_name] = samples_in_group
    
    return expression_df, replicate_groups

def calculate_expression_proportion(gene_expression: pd.Series, threshold: float = 0.5) -> float:
    """
    计算单个基因的表达样本比例。
    """
    expressed_samples = (gene_expression > threshold).sum()
    total_samples = len(gene_expression)
    return expressed_samples / total_samples if total_samples > 0 else 0.0

def calculate_replicate_correlation(gene_expression: pd.Series, replicate_groups: Dict[str, List[str]]) -> float:
    """
    计算单个基因的重复间平均相关性。
    """
    correlations = []
    
    for group_name, samples in replicate_groups.items():
        group_expression = gene_expression[samples].values
        
        if np.all(group_expression == 0):
            continue
            
        n_samples = len(samples)
        if n_samples >= 2:
            for i in range(n_samples):
                for j in range(i + 1, n_samples):
                    corr, _ = pearsonr([group_expression[i]], [group_expression[j]])
                    if not np.isnan(corr):
                        correlations.append(corr)
    
    return np.mean(correlations) if correlations else 0.0

def calculate_replicate_cv(gene_expression: pd.Series, replicate_groups: Dict[str, List[str]]) -> float:
    """
    计算单个基因的重复间平均变异系数(CV)。
    """
    cvs = []
    
    for group_name, samples in replicate_groups.items():
        group_expression = gene_expression[samples].values
        mean_expr = np.mean(group_expression)
        std_expr = np.std(group_expression)
        
        if mean_expr > 0:
            cv = std_expr / mean_expr
            cvs.append(cv)
    
    return np.mean(cvs) if cvs else float('inf')

def calculate_metrics_for_genes(gene_ids: List[str], expression_df: pd.DataFrame, 
                               replicate_groups: Dict[str, List[str]]) -> pd.DataFrame:
    """
    为指定的基因列表计算所有质量指标。
    """
    results = []
    
    for gene_id in gene_ids:
        if gene_id not in expression_df.index:
            continue
            
        gene_expression = expression_df.loc[gene_id]
        
        expr_prop = calculate_expression_proportion(gene_expression)
        replicate_corr = calculate_replicate_correlation(gene_expression, replicate_groups)
        replicate_cv = calculate_replicate_cv(gene_expression, replicate_groups)
        
        results.append({
            'gene_id': gene_id,
            'expression_proportion': expr_prop,
            'replicate_correlation': replicate_corr,
            'replicate_cv': replicate_cv,
            'mean_expression': np.mean(gene_expression)
        })
    
    return pd.DataFrame(results)

def filter_novel_genes(known_gene_ids: List[str], novel_gene_ids: List[str],
                      expression_df: pd.DataFrame, replicate_groups: Dict[str, List[str]],
                      strictness: str = 'medium') -> pd.DataFrame:
    """
    主要筛选函数：基于已知基因建立基准，筛选新基因。
    """
    # 1. 为已知基因计算指标
    print("计算已知基因的质量指标...")
    known_genes_metrics = calculate_metrics_for_genes(known_gene_ids, expression_df, replicate_groups)
    
    # 2. 根据严格度设置分位数
    if strictness == 'lenient':
        quantile = 0.25
    elif strictness == 'strict':
        quantile = 0.1
    else:
        quantile = 0.2
    
    # 3. 建立基准阈值
    thresholds = {
        'min_expression_proportion': known_genes_metrics['expression_proportion'].quantile(quantile),
        'min_replicate_correlation': known_genes_metrics['replicate_correlation'].quantile(quantile),
        'max_replicate_cv': known_genes_metrics['replicate_cv'].quantile(1 - quantile)
    }
    
    print(f"\n建立的基准阈值 ({strictness} 模式):")
    for metric, threshold in thresholds.items():
        print(f"  {metric}: {threshold:.3f}")
    
    # 4. 为新基因计算指标
    print("\n计算新基因的质量指标...")
    novel_genes_metrics = calculate_metrics_for_genes(novel_gene_ids, expression_df, replicate_groups)
    
    # 5. 应用筛选
    mask = (
        (novel_genes_metrics['expression_proportion'] >= thresholds['min_expression_proportion']) &
        (novel_genes_metrics['replicate_correlation'] >= thresholds['min_replicate_correlation']) &
        (novel_genes_metrics['replicate_cv'] <= thresholds['max_replicate_cv'])
    )
    
    novel_genes_metrics['pass_filter'] = mask
    novel_genes_metrics['threshold_expression_proportion'] = thresholds['min_expression_proportion']
    novel_genes_metrics['threshold_replicate_correlation'] = thresholds['min_replicate_correlation']
    novel_genes_metrics['threshold_replicate_cv'] = thresholds['max_replicate_cv']
    
    print(f"\n筛选结果:")
    print(f"  总新基因数: {len(novel_genes_metrics)}")
    print(f"  通过筛选数: {mask.sum()}")
    print(f"  筛选通过率: {mask.sum()/len(novel_genes_metrics)*100:.1f}%")
    
    return novel_genes_metrics

def main():
    """
    主函数：执行完整的筛选流程。
    """
    # 文件路径配置
    expression_matrix_path = "expression_matrix.tsv"  # 修改为您的表达矩阵文件路径
    sample_design_path = "sample_design.tsv"          # 修改为您的实验设计文件路径
    output_path = "novel_genes_expression_screening_results.tsv"
    
    # 加载数据
    print("正在加载数据...")
    expr_df, rep_groups = load_data(expression_matrix_path, sample_design_path)
    print(f"加载完成: {expr_df.shape[0]}个基因, {expr_df.shape[1]}个样本")
    print(f"重复组数: {len(rep_groups)}")
    
    # 获取基因ID列表（这里需要您根据实际情况提供）
    # 已知基因ID列表 - 需要您提供
    known_genes = [gene for gene in expr_df.index if gene.startswith('KNOWN_')]  # 示例，请修改为您的逻辑
    # 新基因ID列表 - 需要您提供  
    novel_genes = [gene for gene in expr_df.index if gene.startswith('NOVEL_')]  # 示例，请修改为您的逻辑
    
    print(f"\n基因统计:")
    print(f"  已知基因数: {len(known_genes)}")
    print(f"  新基因数: {len(novel_genes)}")
    
    if len(known_genes) == 0 or len(novel_genes) == 0:
        print("错误: 未找到足够的已知基因或新基因，请检查基因ID匹配")
        return
    
    # 执行筛选（可选择 'lenient', 'medium', 'strict'）
    results_df = filter_novel_genes(known_genes, novel_genes, expr_df, rep_groups, strictness='medium')
    
    # 保存结果
    results_df.to_csv(output_path, sep='\t', index=False)
    print(f"\n结果已保存至: {output_path}")
    
    # 提取通过筛选的基因
    high_confidence_novel_genes = results_df[results_df['pass_filter']]['gene_id'].tolist()
    print(f"高置信度新基因列表: {len(high_confidence_novel_genes)}个")
    
    # 可选：保存高置信度基因列表
    with open("high_confidence_novel_genes.txt", 'w') as f:
        for gene_id in high_confidence_novel_genes:
            f.write(f"{gene_id}\n")
    print("高置信度新基因ID已保存至: high_confidence_novel_genes.txt")

if __name__ == "__main__":
    main()