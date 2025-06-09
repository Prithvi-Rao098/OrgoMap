'''
import numpy as np
import faiss
from sentence_transformers import SentenceTransformer
from openai import OpenAI
import os

# -------------------------------
# CONFIGURATION
# -------------------------------
INDEX_PATH = "orgo_index.faiss"
DOCS_PATH = "doc_mapping.npy"
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
doc_mapping = np.load(DOCS_PATH, allow_pickle=True)

# -------------------------------
# RAG FUNCTION
# -------------------------------
def ask_rag(query, top_k=TOP_K):
    query_embedding = model.encode([query], convert_to_numpy=True)
    query_embedding = query_embedding / np.linalg.norm(query_embedding)

    D, I = index.search(query_embedding, k=top_k)
    retrieved_context = "\n\n".join([doc_mapping[idx] for idx in I[0]])

    system_prompt = "You are a helpful chemistry tutor. Use the following documents to answer the user's question."
    user_prompt = f"""You are a helpful chemistry tutor. Use the following retrieved documents to answer the question.

Retrieved Documents:
{retrieved_context}

Question: {query}
Answer:"""

    stream = client.chat.completions.create(
        model="gpt-4",
        messages=[
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt}
        ],
        temperature=0.7,
        stream=True,
    )

    print("\n🔬 Answer:")
    for chunk in stream:
        content = chunk.choices[0].delta.content or ""
        print(content, end="", flush=True)

# -------------------------------
# MAIN LOOP
# -------------------------------
if __name__ == "__main__":
    while True:
        try:
            query = input("\n\n📝 Enter your chemistry question: ")
            if query.lower() in ["exit", "quit"]:
                break
            ask_rag(query)
        except KeyboardInterrupt:
            print("\nSession ended.")
            break
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
