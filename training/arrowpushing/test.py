import pandas as pd

file_path = "../../../5104873/1976_Sep2016_USPTOgrants_smiles.rsmi"
sample = pd.read_csv(file_path, sep='\t', header=None, nrows=5)
print("Sample raw data:")
print(sample)