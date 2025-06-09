'''
import pandas as pd
from sentence_transformers import SentenceTransformer
import faiss
import numpy as np
from openai import OpenAI  # New SDK import

print("🚀 Starting script...")

# Initialize OpenAI client with your API key here
client = OpenAI(api_key="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")

# Load CSV
df = pd.read_csv("train.csv")
df.fillna('', inplace=True)
#print(f"✅ CSV Loaded. Number of rows: {len(df)}")

# Combine rows into documents
def create_doc(row):
    return f"""Role: {row['role_1']}
Topic: {row['topic;']}
Sub-topic: {row['sub_topic']}
Question: {row['message_1']}
Answer: {row['message_2']}"""

try:
    docs = df.apply(create_doc, axis=1).tolist()
    #docs = docs[:500]  # For quicker testing
    #print(f"✅ Created documents: {len(docs)}")
    #print("🔍 Sample doc:\n", docs[0])
except Exception as e:
    print(f"❌ Failed creating documents: {e}")
    exit(1)

print("📊 Embedding documents with SentenceTransformer...")
model = SentenceTransformer("all-MiniLM-L6-v2")
embeddings = model.encode(docs, convert_to_numpy=True)
embeddings = embeddings / np.linalg.norm(embeddings, axis=1, keepdims=True)
print("✅ Embeddings created")

print("📚 Creating FAISS index and adding embeddings...")
dimension = embeddings.shape[1]
index = faiss.IndexFlatIP(dimension)
index.add(embeddings)
doc_mapping = {i: doc for i, doc in enumerate(docs)}
print("✅ FAISS index ready")

def ask_rag(query, top_k=3):
    print("🔍 Embedding query...")
    query_embedding = model.encode([query])
    query_embedding = query_embedding / np.linalg.norm(query_embedding)

    print("📚 Searching FAISS index...")
    D, I = index.search(query_embedding, k=top_k)

    retrieved_context = "\n\n".join([doc_mapping[idx] for idx in I[0]])
    print("📄 Retrieved context:")
    #print(retrieved_context)

    prompt = f"""You are a helpful chemistry tutor. Use the following retrieved documents to answer the question.


Retrieved Documents:
{retrieved_context}

Question: {query}
Answer:"""

    print("📨 Sending prompt to GPT-4...")
    try:
        response = client.chat.completions.create(
            model="gpt-4",
            messages=[
                {"role": "system", "content": "You are a helpful chemistry tutor. Use the following documents to help answer the user's question."},
                {"role": "user", "content": prompt}
            ],
            temperature=0.7,
            max_tokens=500,
        )
        answer = response.choices[0].message.content.strip()
        print("✅ GPT-4 Response Received")
        return answer
    except Exception as e:
        print(f"❌ Error calling OpenAI: {e}")
        return "An error occurred while calling the OpenAI API."

if __name__ == "__main__":
    print("🧪 Testing RAG response...")
    answer = ask_rag("How do I name esters?")
    print("🧠 AI Answer:")
    print(answer)
'''

'''
import pandas as pd
import numpy as np
import faiss
from sentence_transformers import SentenceTransformer
from openai import OpenAI
import os
from tqdm import tqdm

# -------------------------------
# CONFIGURATION
# -------------------------------
CSV_PATH = "train.csv"
NUM_DOCS = None  # Set to an integer for limiting rows (e.g., 500), or None for all
BATCH_SIZE = 256
INDEX_PATH = "orgo_index.faiss"
DOCS_PATH = "doc_mapping.npy"

# -------------------------------
# INITIAL SETUP
# -------------------------------
print("🚀 Starting script...")
client = OpenAI(api_key="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")
model = SentenceTransformer("all-MiniLM-L6-v2")

# -------------------------------
# LOAD DATA
# -------------------------------
df = pd.read_csv(CSV_PATH)
df.fillna('', inplace=True)
if NUM_DOCS:
    df = df[:NUM_DOCS]


def create_doc(row):
    return f"""Role: {row['role_1']}
Topic: {row['topic;']}
Sub-topic: {row['sub_topic']}
Question: {row['message_1']}
Answer: {row['message_2']}"""

docs = df.apply(create_doc, axis=1).tolist()

# -------------------------------
# EMBEDDING WITH BATCHING
# -------------------------------
def batched_embeddings(docs, batch_size=BATCH_SIZE):
    embeddings = []
    for i in tqdm(range(0, len(docs), batch_size), desc="Embedding docs"):
        batch = docs[i:i+batch_size]
        batch_embeddings = model.encode(batch, convert_to_numpy=True)
        embeddings.append(batch_embeddings)
    return np.vstack(embeddings)

embeddings = batched_embeddings(docs)
embeddings = embeddings / np.linalg.norm(embeddings, axis=1, keepdims=True)

# -------------------------------
# FAISS INDEX
# -------------------------------
dimension = embeddings.shape[1]
index = faiss.IndexFlatIP(dimension)
index.add(embeddings)

faiss.write_index(index, INDEX_PATH)
with open(DOCS_PATH, "wb") as f:
    np.save(f, docs)

print(f"✅ FAISS index saved to {INDEX_PATH}")
print(f"✅ Documents saved to {DOCS_PATH}")
'''

# data.py
import pandas as pd
import numpy as np
import faiss
from sentence_transformers import SentenceTransformer
from openai import OpenAI
import os
from tqdm import tqdm
import json
import gzip

# -------------------------------
# CONFIGURATION
# -------------------------------
CSV_PATH = "train.csv"
NUM_DOCS = None  # Set to an integer for limiting rows (e.g., 500), or None for all
BATCH_SIZE = 256
INDEX_PATH = "orgo_index.faiss"
DOCS_PATH = "doc_mapping.json.gz"

# -------------------------------
# INITIAL SETUP
# -------------------------------
client = OpenAI(api_key="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")
model = SentenceTransformer("all-MiniLM-L6-v2")

# -------------------------------
# LOAD DATA
# -------------------------------
df = pd.read_csv(CSV_PATH)
df.fillna('', inplace=True)
if NUM_DOCS:
    df = df[:NUM_DOCS]

def create_doc(row):
    return f"""Role: {row['role_1']}
Topic: {row['topic;']}
Sub-topic: {row['sub_topic']}
Question: {row['message_1']}
Answer: {row['message_2']}"""

docs = df.apply(create_doc, axis=1).tolist()

# -------------------------------
# EMBEDDING WITH BATCHING
# -------------------------------
def batched_embeddings(docs, batch_size=BATCH_SIZE):
    embeddings = []
    for i in tqdm(range(0, len(docs), batch_size), desc="Embedding docs"):
        batch = docs[i:i+batch_size]
        batch_embeddings = model.encode(batch, convert_to_numpy=True)
        embeddings.append(batch_embeddings)
    return np.vstack(embeddings)

embeddings = batched_embeddings(docs)
embeddings = embeddings / np.linalg.norm(embeddings, axis=1, keepdims=True)

# -------------------------------
# FAISS INDEX
# -------------------------------
dimension = embeddings.shape[1]
index = faiss.IndexFlatIP(dimension)
index.add(embeddings)

faiss.write_index(index, INDEX_PATH)
with gzip.open(DOCS_PATH, "wt", encoding="utf-8") as f:
    json.dump(docs, f)


