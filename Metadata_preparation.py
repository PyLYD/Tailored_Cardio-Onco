import pandas as pd
import re

file_path = 'C:\CardioOnco\Datasets\Tertiary dataset TCGA BRCA\Metadata\Sample_ID_Metadata_input.xlsx'
df = pd.read_excel(file_path)

# Function to remove last 12 characters using regex
def remove_last_12_chars(sample_id):
    pattern = re.compile(r'^(.*?)-\d{2}[A-Z]-[A-Z0-9]{4}-\d{2}$')
    match = pattern.match(sample_id)
    if match:
        return match.group(1)
    return sample_id

# Function to categorize based on last 3 characters
def categorize_id(sample_id_trimmed):
    if isinstance(sample_id_trimmed, str):
        if sample_id_trimmed.endswith('01A'):
            return 'Patient'
        elif sample_id_trimmed.endswith('11A'):
            return 'Normal'
    return 'Else'

# Apply the function to create new column
df['Sample_ID_Trimmed'] = df['Sample_ID'].apply(remove_last_12_chars)

df['Category'] = df['Sample_ID_Trimmed'].apply(categorize_id)

# Display the DataFrame
#print(df)

#df.to_csv(r'C:\CardioOnco\Datasets\Tertiary dataset TCGA BRCA\Metadata\Metadata_IDs_updated.tsv', sep='\t', index=False)

m1 = r'C:\CardioOnco\Datasets\Tertiary dataset TCGA BRCA\Metadata\Metadata_IDs_updated.tsv'
metadata_to_filter = pd.read_csv(m1, sep="\t", index_col=False)
print(metadata_to_filter)


m2 = r'C:\CardioOnco\Datasets\Tertiary dataset TCGA BRCA\Metadata\Sample_IDs_to_keep.xlsx'
filtered_Sample_ID = pd.read_excel(m2)
print(filtered_Sample_ID)

filtered_metadata = metadata_to_filter.merge(filtered_Sample_ID, on='Sample_ID')
#filtered_metadata = metadata_to_filter[metadata_to_filter['Sample_ID'].isin(filtered_Sample_ID['Sample_ID'])]    #or using "isin"
print(filtered_metadata)

#filtered_metadata.to_csv(r'C:\CardioOnco\Datasets\Tertiary dataset TCGA BRCA\Metadata\Metadata_IDs_PatientVSNormal_final.tsv', sep='\t', index=False)

#Prepare updated metadata for chemo vs non chemo
filtered_metadata_Patient = filtered_metadata[filtered_metadata['Sample_type'] == "Patient"]     #select only patients since there are also normal samples with same short ID (without -01A, -11B)
filtered_metadata_Patient['Patient_ID'] = filtered_metadata_Patient['Patient_ID'].str[:-4]
print(filtered_metadata_Patient)

file_path_2 = 'C:\CardioOnco\Datasets\Tertiary dataset TCGA BRCA\Metadata\Metadata_chemo_vs_no_chemo.tsv'
m3 = pd.read_csv(file_path_2, sep="\t")
print(m3)
m3_chemo_only = m3[m3['Chemo_status'] == "Chemotherapy"]    #We want ot preserve only
print(m3_chemo_only)

df_merged = pd.merge(filtered_metadata_Patient, m3_chemo_only[['Patient_ID', 'Chemo_status']], on='Patient_ID', how='left')
print(df_merged)
df_merged['Chemo_status'].fillna('No_Chemotherapy', inplace=True)
#df_merged.to_csv(r'C:\CardioOnco\Datasets\Tertiary dataset TCGA BRCA\Metadata\Metadata_Chemo_vs_No_chemo_final.tsv', sep='\t', index=False)



#Prepare Metadata for radiotherapy vs non radiotherapy
m4 = r'C:\CardioOnco\Datasets\Tertiary dataset TCGA BRCA\Metadata\Metadata_radiotherapy.xlsx'
metadata_radio = pd.read_excel(m4)
print(metadata_radio)
df_merged_radio = pd.merge(filtered_metadata_Patient, metadata_radio[['Patient_ID', 'TREATMENT_TYPE']], on='Patient_ID', how='left')
df_merged_radio['TREATMENT_TYPE'].fillna('No Radiation Therapy', inplace=True)
print(df_merged_radio)
df_merged_radio = df_merged_radio.drop(columns=['Sample_type'])
df_merged_radio.to_csv(r'C:\CardioOnco\Datasets\Tertiary dataset TCGA BRCA\Metadata\Metadata_Radio_vs_No_Radio_final.tsv', sep='\t', index=False)