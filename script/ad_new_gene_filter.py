import pandas as pd
import os
from collections import defaultdict

def build_group_samples(sample_df):
    group_samples = defaultdict(list)
    cols = {c.lower(): c for c in sample_df.columns}
    sample_col = cols['sample']
    group_col = cols['group']
    for _, row in sample_df.iterrows():
        sample_id = row[sample_col]
        group_name = row[group_col]
        group_samples[group_name].append(sample_id)
    return dict(group_samples)

def prepare_expr_df(expression_df, group_samples):
    cols = [c.lower() for c in expression_df.columns]
    expression_df = expression_df.set_index(expression_df.columns[cols.index('gene')])
    expression_df = expression_df.dropna(axis=1, how='all')
    wanted = set(sum(group_samples.values(), []))
    keep_cols = [c for c in expression_df.columns if str(c) in wanted]
    return expression_df.loc[:, keep_cols]

def calculate_cv_by_group(expression_df, group_samples):
    cv_results = {}
    for group_name, sample_ids in group_samples.items():
        cols = [s for s in sample_ids if s in expression_df.columns]
        group_expr = expression_df.loc[:, cols]
        mean_expr = group_expr.mean(axis=1)
        std_expr = group_expr.std(axis=1) 
        cv = std_expr / (mean_expr + 1e-6)
        cv_results[group_name] = cv
    cv_df = pd.DataFrame(cv_results)
    cv_df.index.name = expression_df.index.name or "gene_id" # 确定行索引
    return cv_df

def selected_new_genes(cv_df, prefix_annot, prefix_new):
    idx = cv_df.index.astype(str)
    mask_annot = idx.str.startswith(prefix_annot)
    mask_new   = idx.str.startswith(prefix_new)
    annot_cv = cv_df.loc[mask_annot]
    new_cv   = cv_df.loc[mask_new]
    thresholds_per_group = annot_cv.max(axis = 0)
    valid_groups = thresholds_per_group.dropna().index.tolist()
    meets_all = new_cv[valid_groups].ge(thresholds_per_group[valid_groups], axis=1).all(axis=1)
    selected = new_cv.loc[meets_all]
    return thresholds_per_group, selected    

def calculate_cv_count_and_filter(expression_file, sample_file, output_dir):
    expression_df = pd.read_excel(expression_file)
    sample_df = pd.read_excel(sample_file)
    # 样本信息与表达量数据的对应处理
    group_samples = build_group_samples(sample_df)
    expression_df = prepare_expr_df(expression_df, group_samples)
    # 计算变异系数
    cv_df = calculate_cv_by_group(expression_df, group_samples)
    # 筛选符合规律的新基因
    thr_per_group, selected_new_all_groups = selected_new_genes(cv_df, prefix_annot = 'evm', prefix_new = 'strg')
    
    os.makedirs(output_dir, exist_ok=True)
    cv_out = os.path.join(output_dir, "gene_cv_by_group.tsv")
    thr_out = os.path.join(output_dir, "annotated_max_cv_thresholds_per_group.tsv")
    sel_out = os.path.join(output_dir, "selected_new_genes_meet_all_groups.tsv")
    cv_df.to_csv(cv_out, sep='\t', index_label='gene_id')
    thr_per_group.to_csv(thr_out, sep='\t', header=['threshold'])
    selected_new_all_groups.to_csv(sel_out, sep='\t', index_label='gene_id')

def main():
    expression_file = ""
    sample_file = ""
    output_dir = ""
    calculate_cv_count_and_filter(expression_file, sample_file, output_dir)

if __name__ == "__main__":
    main()