import pandas as pd
from rdkit import Chem

def load_rsmi(filepath):
    """Load .rsmi file into DataFrame"""
    return pd.read_csv(
        filepath,
        sep='\t',
        header=None,
        names=['ReactionSMILES', 'PatentID', 'Year'],
        usecols=[0, 1, 2] 
    )

df_grants = load_rsmi("1976_Sep2016_USPTOgrants_smiles.rsmi")
df_apps = load_rsmi("2001_Sep2016_USPTOapplications_smiles.rsmi")

df = pd.concat([df_grants, df_apps]).drop_duplicates()
print(f"Total reactions: {len(df)}")

df['IsOrganic'] = df['ReactionSMILES'].apply(
    lambda x: Chem.MolFromSmiles(x.split('>>')[0]) is not None
)
organic_df = df[df['IsOrganic']]
organic_df.to_parquet("uspto_organic.parquet")