from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import pydeseq2
import pandas as pd
import scanpy as sc     #for PCA
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


df = pd.read_csv(r'C:\CardioOnco\Datasets\Tertiary dataset TCGA BRCA\ALL_RSEM_normal_and_patients\ALL_RSEM.tsv', sep="\t", index_col=False, engine='python')
counts = df.set_index('Hugo_Symbol')



"The effective difference between expected counts from RSEM and raw counts from other tools is often minimal at the gene level. "
"Although, RSEM - I believe - is considered more accurate on account of its smarter handling of multimapped reads."
"It's acceptable to convert RSEM's expected counts to a raw count equivalent by rounding to the nearest integer should you need to. "
"Although if you can preserve it as is, then that's always best"

# Transform the count values to integer. Since the UCSC Xena log2-transforms the count data for some unspecified reason.
counts1 = counts.applymap(lambda x: round(((2 ** x) - 1)))
print(counts1)  # Now we have the initial untransformed raw counts and we can proceed to do the DESeq
# Filter out all the counts that are small or equal to 10
counts1 = counts1[counts1.sum(axis = 1) > 10]


# For PyDESeq2, a table transposition is needed. The gene IDs (Ensembl_ID) should be in the columns
transposed_df = counts1.T
print(transposed_df)
print(transposed_df.columns)

# Count table is ready, now we need to integrate the metadata: Healthy control, Tumor...
