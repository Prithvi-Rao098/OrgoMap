'''

from rdkit import Chem
from rdkit.Chem import Draw

def generate_molecule_image(smiles: str, output_path: str = "molecule.png"):
    # Convert SMILES to molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("❌ Invalid SMILES string.")
        return

    # Generate and save image
    img = Draw.MolToImage(mol, size=(300, 300))
    img.save(output_path)
    print(f"✅ Molecule image saved to: {output_path}")

generate_molecule_image("CCO")  # Ethanol

'''

from rdkit import Chem
from rdkit.Chem import AllChem, Draw

mol = Chem.MolFromSmiles("CCO")
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol)
AllChem.UFFOptimizeMolecule(mol)

# Save as image
img = Draw.MolToImage(mol, size=(300, 300))
img.save("molecule.png")
print("✅ Molecule image saved to molecule.png")
