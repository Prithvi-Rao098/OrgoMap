# main.py
from fastapi import FastAPI, Request
from fastapi.responses import JSONResponse
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware
from query import ask_rag
import os

app = FastAPI()

# Mount the static directory (training/static)
static_dir = os.path.join(os.path.dirname(__file__), "static")
app.mount("/static", StaticFiles(directory=static_dir), name="static")

# CORS for frontend
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Replace * with frontend URL in production
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.post("/chat")
async def chat(request: Request):
    body = await request.json()
    query = body.get("query") or body.get("message")
    if not query:
        return JSONResponse(content={"error": "No query provided"}, status_code=400)

    # Updated: Add instruction to avoid unnecessary SMILES generation
    system_prompt = (
        "You are an organic chemistry tutor. "
        "Only generate a SMILES string if the user asks for a molecule or mentions a chemical compound."
    )

    # Call your custom RAG logic with system prompt
    answer, smiles, image_url, sdf_url = ask_rag(query)

    response_data = {
        "response": answer,
    }

    # Only include images if SMILES was generated
    if smiles:
        response_data["smiles"] = smiles
        response_data["image_url"] = image_url
        response_data["sdf_url"] = sdf_url

    return JSONResponse(response_data)
