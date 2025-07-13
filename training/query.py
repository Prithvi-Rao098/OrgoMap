'''
import numpy as np
import faiss
import json
import gzip
import re
import csv
from sentence_transformers import SentenceTransformer
from openai import OpenAI
from rdkit import Chem
from rdkit.Chem import Draw

# -------------------------------
# CONFIG
# -------------------------------
INDEX_PATH = "orgo_index.faiss"
DOCS_PATH = "doc_mapping.json.gz"
MECH_CSV = "mechanism.csv"
TOP_K = 3

# -------------------------------
# SETUP
# -------------------------------
client = OpenAI(api_key="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")

model = SentenceTransformer("all-MiniLM-L6-v2")
index = faiss.read_index(INDEX_PATH)

with gzip.open(DOCS_PATH, "rt", encoding="utf-8") as f:
    docs = json.load(f)

# -------------------------------
# LOAD STATIC IMAGE MAP
# -------------------------------
static_image_map = {}
with open(MECH_CSV, newline='', encoding='utf-8') as csvfile:
    reader = csv.DictReader(csvfile)
    fieldnames = reader.fieldnames
    print("CSV fields:", fieldnames)
    for row in reader:
        smiles = row['SMILES'].strip()
        path = row[fieldnames[0]].strip()  # first column is file path
        static_image_map[smiles] = path

print(f"✅ Loaded {len(static_image_map)} static images")

# -------------------------------
# DRAW SMILES IMAGE
# -------------------------------
def draw_smiles_image(smiles, filename="molecule.png"):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        img = Draw.MolToImage(mol, size=(300, 300))
        img.save(filename)
        print(f"🧪 2D structure saved as '{filename}'")
        return True
    return False

# -------------------------------
# EXTRACT SMILES
# -------------------------------
def extract_smiles(text):
    match = re.search(r"SMILES\s*:\s*([A-Za-z0-9@+\-\[\]\(\)=#$*]+)", text)
    if match:
        return match.group(1)
    return None

# -------------------------------
# ASK + DRAW LOGIC
# -------------------------------
def ask_rag(query):
    query_embedding = model.encode([query], convert_to_numpy=True)
    query_embedding /= np.linalg.norm(query_embedding, axis=1, keepdims=True)
    D, I = index.search(query_embedding, TOP_K)
    context = "\n\n".join([docs[i] for i in I[0]])

    messages = [
        {
            "role": "system",
            "content": (
                "You are a precise chemistry assistant.\n"
                "- If a molecule is described, always output the *exact SMILES* you know for it.\n"
                "- If the SMILES is unusual or incomplete, still output it exactly.\n"
                "- Always end your answer with: SMILES: <smiles_string>\n"
                "- If no structure applies, end with: SMILES: None"
            )
        },
        {
            "role": "user",
            "content": f"Use the following context to answer the question.\n\nContext:\n{context}\n\nQuestion: {query}"
        }
    ]

    print("\n🔬 Answer:")
    full_answer = ""
    response = client.chat.completions.create(
        model="gpt-4o",
        messages=messages,
        stream=True,
    )
    for chunk in response:
        delta = chunk.choices[0].delta
        if hasattr(delta, "content") and delta.content:
            print(delta.content, end="", flush=True)
            full_answer += delta.content

    print("\n---------------------")
    smiles = extract_smiles(full_answer)
    if smiles and smiles != "None":
        print(f"✅ Found SMILES: {smiles}")
        if "*" in smiles or "[" in smiles:
            # Wildcard or atom-mapping -> try static fallback
            static_path = static_image_map.get(smiles)
            if static_path:
                print(f"📷 Using static image: {static_path}")
            else:
                print(f"⚠️ Wildcard SMILES but no static image found for: {smiles}")
        else:
            if not draw_smiles_image(smiles):
                print("⚠️ RDKit could not draw it. No image generated.")
    else:
        print("⚠️ No SMILES found. No image generated.")

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
            print("\n👋 Session ended.")
            break
'''

import os
import numpy as np
import faiss
import json
import gzip
import re
import csv
from sentence_transformers import SentenceTransformer
from openai import OpenAI
from rdkit import Chem
from rdkit.Chem import Draw, AllChem

# -------------------------------
# CONFIG
# -------------------------------
INDEX_PATH = "orgo_index.faiss"
DOCS_PATH = "doc_mapping.json.gz"
MECH_CSV = "mechanism.csv"
TOP_K = 3

OUTPUT_2D = "static/molecule.png"
OUTPUT_3D = "static/molecule.mol"

os.makedirs("static", exist_ok=True)

# -------------------------------
# SETUP
# -------------------------------
client = OpenAI(api_key="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")
model = SentenceTransformer("all-MiniLM-L6-v2")
index = faiss.read_index(INDEX_PATH)

with gzip.open(DOCS_PATH, "rt", encoding="utf-8") as f:
    docs = json.load(f)

# -------------------------------
# STATIC IMAGE MAP
# -------------------------------
static_image_map = {}
with open(MECH_CSV, newline='', encoding='utf-8') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        smiles = row['SMILES'].strip()
        path = row[reader.fieldnames[0]].strip()
        static_image_map[smiles] = path
print(f"✅ Loaded {len(static_image_map)} static images")

# -------------------------------
# Helpers
# -------------------------------
def draw_smiles_image(smiles, filename=OUTPUT_2D):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        img = Draw.MolToImage(mol, size=(300, 300))
        img.save(filename)
        print(f"🧪 2D structure saved: {filename}")
        return True
    return False

def generate_3d_structure(smiles, filename=OUTPUT_3D):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        mol = Chem.AddHs(mol)
        success = AllChem.EmbedMolecule(mol)
        if success == 0:
            AllChem.UFFOptimizeMolecule(mol)
            Chem.MolToMolFile(mol, filename)
            print(f"🔷 3D structure saved: {filename}")
            return True
    return False

def extract_smiles(text):
    match = re.search(r"SMILES\s*:\s*([A-Za-z0-9@+\-\[\]\(\)=#$*]+)", text)
    if match:
        return match.group(1)
    return None

# -------------------------------
# Main logic
# -------------------------------
def ask_rag(query):
    query_embedding = model.encode([query], convert_to_numpy=True)
    query_embedding /= np.linalg.norm(query_embedding, axis=1, keepdims=True)
    D, I = index.search(query_embedding, TOP_K)
    context = "\n\n".join([docs[i] for i in I[0]])

    messages = [
        {
            "role": "system",
            "content": (
                "You are a precise chemistry assistant.\n"
                "- Always output the exact SMILES.\n"
                "- Always end with: SMILES: <smiles_string> or SMILES: None"
            )
        },
        {
            "role": "user",
            "content": f"Context:\n{context}\n\nQuestion: {query}"
        }
    ]

    print("\n🔬 Answer:")
    full_answer = ""
    response = client.chat.completions.create(
        model="gpt-4o",
        messages=messages,
        stream=True,
    )
    for chunk in response:
        delta = chunk.choices[0].delta
        if hasattr(delta, "content") and delta.content:
            print(delta.content, end="", flush=True)
            full_answer += delta.content

    print("\n---------------------")
    smiles = extract_smiles(full_answer)
    if smiles and smiles != "None":
        print(f"✅ Found SMILES: {smiles}")

        if "*" in smiles or "[" in smiles:
            static_path = static_image_map.get(smiles)
            if static_path:
                print(f"📷 Using static image: {static_path}")
            else:
                print(f"⚠️ Wildcard SMILES but no static image found.")
        else:
            drew = draw_smiles_image(smiles)
            made3d = generate_3d_structure(smiles)
            if not drew:
                print("⚠️ RDKit could not draw 2D.")
            if not made3d:
                print("⚠️ RDKit could not generate 3D.")
            print(f"👉 Serve `{OUTPUT_2D}` for 2D or `{OUTPUT_3D}` for 3D.")
    else:
        print("⚠️ No SMILES found. Nothing to render.")

# -------------------------------
# MAIN LOOP
# -------------------------------
if __name__ == "__main__":
    while True:
        try:
            query = input("\n📝 Enter your chemistry question (or 'exit' to quit): ")
            if query.lower() in ["exit", "quit"]:
                print("👋 Goodbye!")
                break
            ask_rag(query)
        except KeyboardInterrupt:
            print("\n👋 Session ended.")
            break

