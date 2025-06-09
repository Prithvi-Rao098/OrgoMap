import pandas as pd
from sentence_transformers import SentenceTransformer
import faiss
import numpy as np

# Load CSV
df = pd.read_csv("train.csv")
df.fillna('', inplace=True)

# Build text entries for embedding
def create_doc(row):
    return f"""Role: {row['role_1']}
Topic: {row['topic']}
Sub-topic: {row['sub_topic']}
Question: {row['message']}
Answer: {row['message_2']}"""

docs = df.apply(create_doc, axis=1).tolist()

# Embed the documents
model = SentenceTransformer('all-MiniLM-L6-v2')
embeddings = model.encode(docs, convert_to_numpy=True)
embeddings = embeddings / np.linalg.norm(embeddings, axis=1, keepdims=True)

# Create FAISS index
dimension = embeddings.shape[1]
index = faiss.IndexFlatIP(dimension)
index.add(embeddings)

# Optional: Store a mapping to original content
doc_mapping = {i: doc for i, doc in enumerate(docs)}

# Sample query
query = "How do I name esters?"
query_embedding = model.encode([query])
query_embedding = query_embedding / np.linalg.norm(query_embedding)

# Search for the top 3 most relevant docs
D, I = index.search(query_embedding, k=3)
for idx in I[0]:
    print(doc_mapping[idx])
