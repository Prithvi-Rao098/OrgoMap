import os
from generation import make_diagram  # Assuming your main code is in reaction_diagram.py

def generate_test_questions():
    # Create test_questions directory if it doesn't exist
    os.makedirs("test_questions", exist_ok=True)
    
    # List of test questions (reaction SMILES)
    test_reactions = [
        # Basic reactions
        "CCO>>C=O",                          # Alcohol oxidation
        "C=O.[H][H]>>CO",                    # Carbonyl reduction
        "C=C.C=C>>C1CCCC1",                  # Diels-Alder
        "CBr>>CO",                           # Substitution
        "CC(=O)OC>>CC(=O)O",                 # Ester hydrolysis
        
        # More complex examples
        "[Cl:1][C:2]([H:3])([H:4])[H:5].[OH2:6]>>[O:6][C:2]([H:3])([H:4])[H:5].[Cl:1][H:6]",  # Hydrolysis
        "[C:1](=[O:2])[OH:3]>>[C:1](=[O:2])[O:3][C:4]",  # Esterification
        "[N:1]([H:2])([H:3])[C:4]>>[N+:1]([H:2])([H:3])[C:4]",  # Amine protonation
        
        # Named reactions
        "C1=CC=CC=C1.C=O>>C1=CC=CC=C1C=O",  # Friedel-Crafts acylation
        "C#N>>C=O",                         # Nitrile hydrolysis
    ]
    
    # Generate diagrams for each test question
    for i, reaction_smiles in enumerate(test_reactions, 1):
        try:
            output_file = f"test_questions/reaction_{i}.svg"
            make_diagram(reaction_smiles, out_svg=output_file, show_arrows=False)
            print(f"Generated: {output_file}")
        except Exception as e:
            print(f"Failed to generate reaction {i}: {e}")

if __name__ == "__main__":
    generate_test_questions()
    print("All test question diagrams generated in 'test_questions' directory")