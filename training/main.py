# main.py
from fastapi import FastAPI, Request
from fastapi.responses import JSONResponse, FileResponse
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware
from query import ask_rag
import os
import time



app = FastAPI()

# This resolves to training/static/
static_dir = os.path.join(os.path.dirname(__file__), "static")
app.mount("/static", StaticFiles(directory=static_dir), name="static")


# Allow frontend
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # You can restrict this to your frontend origin
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Serve static files (molecule.png, molecule.sdf)
app.mount("/static", StaticFiles(directory="static"), name="static")

@app.post("/chat")  # change from /ask_rag to /chat
async def chat(request: Request):
    body = await request.json()
    query = body.get("query") or body.get("message")  # handle both keys if needed
    if not query:
        return JSONResponse(content={"error": "No query provided"}, status_code=400)

    answer, smiles = ask_rag(query)
    return JSONResponse({
        "response": answer,
        "smiles": smiles,
        "image_url": f"http://localhost:8000/static/molecule.png?t={int(time.time())}",
        "sdf_url": "http://localhost:8000/static/molecule.sdf"
    })

