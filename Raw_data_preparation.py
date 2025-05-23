import pandas as pd

ALL_RSEM = r'C:\CardioOnco\Datasets\Tertiary dataset TCGA BRCA\ALL_RSEM_normal_and_patients\BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.tsv'
df10 = pd.read_csv(ALL_RSEM, sep="\t", index_col=False)      #at the beginning, take only a sample for testing nrows = 3
#print(df1)
Metadata = r"C:\CardioOnco\Datasets\Tertiary dataset TCGA BRCA\Metadata\Metadata_PatientVSnormal.tsv"
df2 = pd.read_csv(Metadata, sep="\t")
#print(df2.columns)
#Convert sample ID with Patient ID

#print(df1)

df_transposed = df10.T
#print(df_transposed)

df_filtered = df_transposed[(df_transposed.iloc[:, 0] == "raw_count") | (df_transposed.iloc[:, 0] == "gene_id")]
#print(df_filtered)

df_transposed_2 = df_filtered.T
#print(df_transposed_2)

#STOPPED HERE, EXPORT THE DF TO TSV IN ORDER TO RE RUN THE DEA EVENTUALLY, NOW WE ARE SURE WE ARE WORKING WITH THE RAW COUNTS (check also GDC if the count values are the same)...
#Counts in cBioPortal were log transformed based on the 75th percentile of the raw data. It was impossible to retrieve the raw counts from that dataset and it was useless because DESEQ2 accepts only raw counts
#NEXT IS TO REPLACE THE GENE IDs with the readable HUGO symbols

#df_transposed_2.to_csv(r'C:\CardioOnco\Datasets\Tertiary dataset TCGA BRCA\ALL_RSEM_normal_and_patients\BRCA_ALL_RSEM_raw_counts_v3.tsv', sep='\t', index=False)


ALL_RSEM_2 = r'C:\CardioOnco\Datasets\Tertiary dataset TCGA BRCA\ALL_RSEM_normal_and_patients\BRCA_ALL_RSEM_raw_counts_v3.tsv'
df20 = pd.read_csv(ALL_RSEM_2, sep="\t", index_col=False)
df20['Hugo_Symbol'] = df20['gene_id'].str.split('|').str[0]
print(df20)

df20.to_csv(r'C:\CardioOnco\Datasets\Tertiary dataset TCGA BRCA\ALL_RSEM_normal_and_patients\BRCA_ALL_RSEM_raw_counts_final.tsv', sep='\t', index=False)