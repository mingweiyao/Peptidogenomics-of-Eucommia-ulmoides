import pandas as pd
input_file = "/Users/lemon/Desktop/output_with_no_predict/filter_by_ms_transcripts_with_peptides.gtf"
df = pd.read_csv(
        input_file, sep='\t', comment='#', header=None,
        names=['chrom','source','feature','start','end','score','strand','frame','attrs'],
        dtype={'chrom':'string','source':'string','feature':'string','start':'int64','end':'int64','score':'string','strand':'string','frame':'string','attrs':'string'},
        engine='c'
    )
feature_type = df[df['feature'] == 'transcript']
print(len(feature_type))
