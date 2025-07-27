import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor
from rdkit.Chem import rdChemReactions as Reactions
from rdkit.Chem.Draw import rdMolDraw2D as Draw
from rdkit.Geometry import Point2D
import math
import os

HETERO = {7, 8, 15, 16, 17, 35, 53}  # N,O,P,S,Cl,Br,I
HALO = {9, 17, 35, 53}               # F,Cl,Br,I

# -------- parsing --------
def rxn_from_smiles(rsmi: str):
    """Parse reaction SMILES robustly."""
    rxn = Reactions.ReactionFromSmarts(rsmi, useSmiles=True)
    if rxn is None and '>' in rsmi:           # drop agents if needed
        parts = rsmi.split('>')
        if len(parts) == 3:
            rxn = Reactions.ReactionFromSmarts(parts[0] + '>>' + parts[2], useSmiles=True)
    return rxn

def split_rxn(rxn):
    R = [Chem.Mol(rxn.GetReactantTemplate(i)) for i in range(rxn.GetNumReactantTemplates())]
    P = [Chem.Mol(rxn.GetProductTemplate(i))  for i in range(rxn.GetNumProductTemplates())]
    return R, P

# -------- bond-change detection --------
def bonds_by_map(mol):
    s = set()
    for b in mol.GetBonds():
        a, c = b.GetBeginAtom(), b.GetEndAtom()
        ma, mc = a.GetAtomMapNum(), c.GetAtomMapNum()
        if ma and mc:
            s.add(tuple(sorted((ma, mc))))
    return s

def union_bonds(mols):
    out = set()
    for m in mols:
        out |= bonds_by_map(m)
    return out

def has_mapped_to_unmapped_halogen(reactants):
    for m in reactants:
        for b in m.GetBonds():
            a, c = b.GetBeginAtom(), b.GetEndAtom()
            ma, mc = a.GetAtomMapNum(), c.GetAtomMapNum()
            Za, Zc = a.GetAtomicNum(), c.GetAtomicNum()
            if (ma and not mc and Zc in HALO) or (mc and not ma and Za in HALO):
                return True
    return False

# -------- drawing functions --------
def draw_clean_reaction(reactants, products, fn="clean_reaction.svg"):
    """Draws reaction without arrows or atom numbers"""
    # Create copies to avoid modifying original molecules
    reactants = [Chem.Mol(mol) for mol in reactants]
    products = [Chem.Mol(mol) for mol in products]
    
    # Remove atom map numbers
    for mol in reactants + products:
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
    
    # Create reaction object
    rxn = AllChem.ChemicalReaction()
    for mol in reactants: rxn.AddReactantTemplate(mol)
    for mol in products: rxn.AddProductTemplate(mol)
    
    # Set up drawer
    drawer = Draw.MolDraw2DSVG(800, 400)
    dopts = drawer.drawOptions()
    
    # Customize appearance
    dopts.addAtomIndices = False  # Hide atom numbers
    dopts.includeAtomTags = False # Hide additional labels
    dopts.bondLineWidth = 2       # Thicker bonds
    
    # Draw reaction
    drawer.DrawReaction(rxn)
    
    drawer.FinishDrawing()
    with open(fn, 'w') as f:
        f.write(drawer.GetDrawingText())
    return fn

def draw_reaction_with_arrows(reactants, products, fn="reaction_with_arrows.svg"):
    """Draws reaction with electron-pushing arrows"""
    # First draw the clean reaction
    clean_fn = draw_clean_reaction(reactants, products, "temp.svg")
    
    # Then add arrows (using your existing arrow drawing code)
    # ... [your arrow drawing implementation here] ...
    
    # For now we'll just return the clean version
    os.rename("temp.svg", fn)
    return fn

# -------- API --------
def make_diagram(reaction_smiles, out_svg="reaction.svg", show_arrows=False):
    """Main function to generate reaction diagrams"""
    rxn = rxn_from_smiles(reaction_smiles)
    if not rxn:
        raise ValueError("Could not parse reaction SMILES.")
    
    R, P = split_rxn(rxn)
    
    if show_arrows:
        return draw_reaction_with_arrows(R, P, out_svg)
    else:
        return draw_clean_reaction(R, P, out_svg)

def batch_generate_diagrams(parquet_path, output_dir, max_reactions=1000, show_arrows=False):
    """Generate diagrams for multiple reactions from a parquet file"""
    os.makedirs(output_dir, exist_ok=True)
    
    df = pd.read_parquet(parquet_path)
    successful = 0
    
    for i, row in df.head(max_reactions).iterrows():
        try:
            rsmi = row['ReactionSMILES_canonical'] or row['ReactionSMILES']
            out_path = os.path.join(output_dir, f"reaction_{i}.svg")
            make_diagram(rsmi, out_path, show_arrows)
            successful += 1
        except Exception as e:
            print(f"Failed to process row {i}: {e}")
    
    print(f"Successfully generated {successful}/{min(max_reactions, len(df))} diagrams")

# -------- main script --------
if __name__ == "__main__":
    # Example usage for batch processing without arrows
    batch_generate_diagrams(
        "./uspto_parquet/1976_Sep2016_USPTOgrants_smiles.rsmi.part0.parquet",
        "./reaction_diagrams",
        max_reactions=1000,
        show_arrows=False  # Set to True if you want arrows
    )