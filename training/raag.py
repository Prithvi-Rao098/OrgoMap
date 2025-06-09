from vertexai.preview.generative_models import GenerativeModel
import faiss
import pickle

# Load index and answer mapping
index = faiss.read_index("chem_index.faiss")
with open("chem_answers.pkl", "rb") as f:
    answers = pickle.load(f)

# Load embedder again
embedder = SentenceTransformer("all-MiniLM-L6-v2")

# Initialize Gemini 1.5 Pro
model = GenerativeModel("gemini-1.5-pro")

def ask_chem_question(user_question, top_k=3):
    # Embed the user's question
    query_embedding = embedder.encode([user_question])
    
    # Retrieve top matches
    D, I = index.search(np.array(query_embedding), top_k)
    context = "\n".join([f"Q: {questions[i]}\nA: {answers[i]}" for i in I[0]])

    # Construct Gemini prompt
    prompt = f"""You are a chemistry tutor. Use the following context to answer the question.
{context}

User Question: {user_question}
Answer:"""

    response = model.generate_content(prompt)
    return response.text
