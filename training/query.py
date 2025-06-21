'''

# query.py
import faiss
import numpy as np
import json
import gzip
from sentence_transformers import SentenceTransformer
from openai import OpenAI

# -------------------------------
# CONFIGURATION
# -------------------------------
INDEX_PATH = "orgo_index.faiss"
DOCS_PATH = "doc_mapping.json.gz"
TOP_K = 3

# -------------------------------
# INITIAL SETUP
# -------------------------------
client = OpenAI(api_key="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")
model = SentenceTransformer("all-MiniLM-L6-v2")

# -------------------------------
# LOAD INDEX AND DOCUMENTS
# -------------------------------
index = faiss.read_index(INDEX_PATH)

with gzip.open(DOCS_PATH, "rt", encoding="utf-8") as f:
    docs = json.load(f)

# -------------------------------
# ASK FUNCTION
# -------------------------------
def ask_rag(query):
    query_embedding = model.encode([query], convert_to_numpy=True)
    query_embedding /= np.linalg.norm(query_embedding, axis=1, keepdims=True)

    D, I = index.search(query_embedding, TOP_K)
    context = "\n\n".join([docs[i] for i in I[0]])

    messages = [
        {"role": "system", "content": "You are a helpful chemistry assistant."},
        {"role": "user", "content": f"Use the following context to answer the question.\n\nContext:\n{context}\n\nQuestion: {query}"}
    ]

    print("\n\n🔬 Answer:")
    response = client.chat.completions.create(
        model="gpt-4",
        messages=messages,
        stream=True,
    )
    for chunk in response:
        delta = chunk.choices[0].delta
        if hasattr(delta, 'content') and delta.content:
            print(delta.content, end='', flush=True)

# -------------------------------
# MAIN LOOP
# -------------------------------
if __name__ == "__main__":
    query = input("\n\n📝 Enter your chemistry question: ")
    ask_rag(query)
    
'''
import numpy as np
import faiss
import json
import gzip
import re
from sentence_transformers import SentenceTransformer
from openai import OpenAI
from rdkit import Chem
from rdkit.Chem import Draw

# -------------------------------
# CONFIGURATION
# -------------------------------
INDEX_PATH = "orgo_index.faiss"
DOCS_PATH = "doc_mapping.json.gz"
TOP_K = 3

# -------------------------------
# INITIAL SETUP
# -------------------------------
client = OpenAI(api_key="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")
model = SentenceTransformer("all-MiniLM-L6-v2")

# -------------------------------
# LOAD INDEX AND DOCUMENTS
# -------------------------------
index = faiss.read_index(INDEX_PATH)

with gzip.open(DOCS_PATH, "rt", encoding="utf-8") as f:
    docs = json.load(f)

# -------------------------------
# HELPER: Draw and save molecule image
# -------------------------------
def draw_smiles_image(smiles, filename="molecule.png"):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        img = Draw.MolToImage(mol, size=(300, 300))
        img.save(filename)
        print(f"\n🧪 2D structure saved as '{filename}'")
    else:
        print("\n⚠️ Could not parse SMILES string for drawing.")

# -------------------------------
# EXTRACT SMILES FROM TEXT
# -------------------------------
def extract_smiles(text):
    # Look for a line like "SMILES: CCO"
    match = re.search(r"SMILES\s*:\s*([A-Za-z0-9@+\-\[\]\(\)=#$]+)", text)
    if match:
        return match.group(1)
    return None

# -------------------------------
# ASK FUNCTION WITH STREAMING AND IMAGE GENERATION
# -------------------------------
def ask_rag(query):
    query_embedding = model.encode([query], convert_to_numpy=True)
    query_embedding /= np.linalg.norm(query_embedding, axis=1, keepdims=True)

    D, I = index.search(query_embedding, TOP_K)
    context = "\n\n".join([docs[i] for i in I[0]])

    messages = [
        {"role": "system", "content": (
            "You are a helpful chemistry assistant. "
            "If you mention a chemical structure, always provide the SMILES notation "
            "in the format: SMILES: <smiles_string>."
        )},
        {"role": "user", "content": f"Use the following context to answer the question.\n\nContext:\n{context}\n\nQuestion: {query}"}
    ]

    print("\n\n🔬 Answer:")
    full_answer = ""
    response = client.chat.completions.create(
        model="gpt-4",
        messages=messages,
        stream=True,
    )
    for chunk in response:
        delta = chunk.choices[0].delta
        if hasattr(delta, 'content') and delta.content:
            print(delta.content, end='', flush=True)
            full_answer += delta.content

    # After full answer is printed, extract SMILES and draw image if found
    smiles = extract_smiles(full_answer)
    if smiles:
        draw_smiles_image(smiles)
    else:
        print("\n⚠️ No valid SMILES notation found in the answer.")

# -------------------------------
# MAIN LOOP
# -------------------------------
if __name__ == "__main__":
    while True:
        try:
            query = input("\n\n📝 Enter your chemistry question (or 'exit' to quit): ")
            if query.lower() in ["exit", "quit"]:
                print("👋 Goodbye!")
                break
            ask_rag(query)
        except KeyboardInterrupt:
            print("\nSession ended.")
            break

