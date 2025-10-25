import pandas as pd
import os

def _build_group_samples(metadata_df: pd.DataFrame) -> dict:
    cols = {c.lower(): c for c in metadata_df.columns}
    sample_col = cols['sample']; group_col = cols['group']
    group_samples = {}
    for _, row in metadata_df.iterrows():
        sample_id = str(row[sample_col])
        group_name = str(row[group_col])
        group_samples.setdefault(group_name, []).append(sample_id)
    return group_samples

def _prepare_expr_df(expr_df: pd.DataFrame, group_samples: dict) -> pd.DataFrame:
    lower_cols = [c.lower() for c in expr_df.columns]
    if 'gene' in lower_cols:
        expr_df = expr_df.set_index(expr_df.columns[lower_cols.index('gene')])
    elif 'gene_id' in lower_cols:
        expr_df = expr_df.set_index(expr_df.columns[lower_cols.index('gene_id')])
    expr_df = expr_df.dropna(axis=1, how='all')
    wanted = set(sum(group_samples.values(), []))
    keep_cols = [c for c in expr_df.columns if str(c) in wanted]
    if len(keep_cols) == 0:
        raise ValueError("表达矩阵里没有与元数据匹配的样本列。请检查样本ID是否一致。")
    return expr_df.loc[:, keep_cols]

def _calculate_cv_by_group(expr_df: pd.DataFrame, group_samples: dict) -> pd.DataFrame:
    cv_results = {}
    for group_name, sample_ids in group_samples.items():
        cols = [s for s in sample_ids if s in expr_df.columns]
        if len(cols) == 0:
            cv_results[group_name] = pd.Series(float('nan'), index=expr_df.index)
            continue
        group_expr = expr_df.loc[:, cols]
        mean_expr = group_expr.mean(axis=1)
        std_expr = group_expr.std(axis=1) 
        cv = std_expr / (mean_expr + 1e-6)
        cv_results[group_name] = cv
    cv_df = pd.DataFrame(cv_results)
    cv_df.index.name = expr_df.index.name or "gene_id"
    return cv_df

def _select_new_genes_meet_all_group_thresholds(cv_df: pd.DataFrame,
                                                prefix_annot='evm',
                                                prefix_new='strg'):
    idx = cv_df.index.astype(str)
    mask_annot = idx.str.startswith(prefix_annot)
    mask_new   = idx.str.startswith(prefix_new)
    annot_cv = cv_df.loc[mask_annot]
    new_cv   = cv_df.loc[mask_new]
    thresholds_per_group = annot_cv.max(axis=0)
    valid_groups = thresholds_per_group.dropna().index.tolist()
    if len(valid_groups) == 0:
        selected = new_cv.iloc[0:0].copy()
        return thresholds_per_group, selected
    meets_all = new_cv[valid_groups].ge(thresholds_per_group[valid_groups], axis=1).all(axis=1)
    selected = new_cv.loc[meets_all].sort_values(by=valid_groups, ascending=False)
    return thresholds_per_group, selected

def calculate_cv_count_by_group(expression_file, metadata_file, output_dir):
    expr_df = pd.read_excel(expression_file)
    metadata_df = pd.read_excel(metadata_file)

    group_samples = _build_group_samples(metadata_df)
    expr_df = _prepare_expr_df(expr_df, group_samples)
    cv_df = _calculate_cv_by_group(expr_df, group_samples)
    thr_per_group, selected_new_all_groups = _select_new_genes_meet_all_group_thresholds(
        cv_df, prefix_annot='evm', prefix_new='strg'
    )

    os.makedirs(output_dir, exist_ok=True)
    cv_out = os.path.join(output_dir, "gene_cv_by_group.tsv")
    thr_out = os.path.join(output_dir, "annotated_max_cv_thresholds_per_group.tsv")
    sel_out = os.path.join(output_dir, "selected_new_genes_meet_all_groups.tsv")
    cv_df.to_csv(cv_out, sep='\t', index_label='gene_id')
    thr_per_group.to_csv(thr_out, sep='\t', header=['threshold'])
    selected_new_all_groups.to_csv(sel_out, sep='\t', index_label='gene_id')

def main():
    expression_file = "path/to/your/expression_matrix.xlsx"
    metadata_file   = "path/to/your/sample_metadata.xlsx"
    output_dir      = "cv_analysis_results"
    calculate_cv_count_by_group(expression_file, metadata_file, output_dir)

if __name__ == "__main__":
    main()
