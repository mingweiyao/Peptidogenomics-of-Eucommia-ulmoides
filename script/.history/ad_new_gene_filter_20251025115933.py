import pandas as pd
import os

def calculate_cv_count_by_group(expression_file, metadata_file, output_dir):
    expr_df = pd.read_excel(expression_file)
    metadata_df = pd.read_excel(metadata_file)
    group_samples = {}

    for _, row in metadata_df.iterrows():
        sample_id = row['sample']
        group_name = row['group']
        if group_name not in group_samples:
            group_samples[group_name] = []
        group_samples[group_name].append(sample_id)

    cv_results = {}
    for group_name, sample_ids in group_samples.items():
        group_expr = expr_df[sample_ids]
        mean_expr = group_expr.mean(axis=1)
        std_expr = group_expr.std(axis=1)
        cv = std_expr / (mean_expr + 1e-6)
        cv_results[group_name] = cv

    cv_df = pd.DataFrame(cv_results)
    cv_df.index = expr_df.index

    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, "gene_cv_by_group.tsv")
    cv_df.to_csv(output_file, sep='\t')
    print(f"CV结果已保存至: {output_file}")

def main():
    expression_file = "path/to/your/expression_matrix.xlsx"
    metadata_file = "path/to/your/sample_metadata.xlsx"
    output_dir = "cv_analysis_results"
    calculate_cv_count_by_group(expression_file, metadata_file, output_dir)

if __name__ == "__main__":
    main()