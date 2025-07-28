import os
import re
from openai import OpenAI
from rdkit import Chem
from dotenv import load_dotenv
from generation import make_diagram

load_dotenv()
client = OpenAI(api_key=os.getenv("VITE_OPENAI_API_KEY"))

# Predefined reactions to ensure validity
PREDEFINED_REACTIONS = {
    "aldol condensation between acetaldehyde and benzaldehyde": "O=CC.[CH]=O.c1ccccc1>>O=CC=CC(=O)c1ccccc1",
    "reaction of hbr with propene": "C=CC.Br>>CC(Br)C",
    "electrophilic aromatic substitution of bromobenzene": "Brc1ccccc1.O=N(=O)O>>Brc1ccc([N+](=O)[O-])cc1"
}

def validate_smiles(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        Chem.SanitizeMol(mol)
        return True
    except:
        return False

def query_to_smiles(query):
    # Check if query matches any predefined reaction
    query_lower = query.lower()
    for pattern, smiles in PREDEFINED_REACTIONS.items():
        if pattern in query_lower:
            return smiles
    
    prompt = f"""Generate valid reaction SMILES for: {query}

Rules:
1. Use standard SMILES syntax
2. For aromatic rings, use lowercase 'c' not 'C'
3. Properly represent charges (e.g., [N+] not [N+])
4. Format as: Reactants>>Products
5. Example: "C=CC.Br>>CC(Br)C" for HBr + propene
"""

    response = client.chat.completions.create(
        model="gpt-4",
        messages=[
            {"role": "system", "content": "Generate valid reaction SMILES following all rules."},
            {"role": "user", "content": prompt}
        ],
        temperature=0.3
    )
    
    smiles = response.choices[0].message.content.strip()
    if not validate_reaction_smiles(smiles):
        raise ValueError(f"Generated invalid SMILES: {smiles}")
    return smiles

def validate_reaction_smiles(smiles):
    if '>>' not in smiles:
        return False
    
    reactants, products = smiles.split('>>')
    for part in [reactants, products]:
        for component in part.split('.'):
            if not validate_smiles(component):
                return False
    return True

def generate_diagram(query):
    os.makedirs("reaction_images", exist_ok=True)
    
    try:
        print(f"\nProcessing: {query}")
        
        reaction_smiles = query_to_smiles(query)
        print(f"Using SMILES: {reaction_smiles}")
        
        # Fixed filename generation - moved regex pattern outside f-string
        clean_query = re.sub(r'[^\w]', '_', query.lower())[:50]
        filename = f"reaction_images/{clean_query}.png"
        
        temp_svg = "temp_rxn.svg"
        make_diagram(reaction_smiles, temp_svg)
        
        import cairosvg
        cairosvg.svg2png(url=temp_svg, write_to=filename)
        os.remove(temp_svg)
        
        print(f"Saved to: {filename}")
        return filename
        
    except Exception as e:
        print(f"Failed: {str(e)}")
        return None

if __name__ == "__main__":
    print("Reaction Diagram Generator (type 'quit' to exit)")
    while True:
        query = input("\nDescribe reaction: ").strip()
        if query.lower() in ('quit', 'exit', 'q'):
            break
        if query:
            generate_diagram(query)