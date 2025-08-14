# server/main.py

import os
import re
import hashlib
from typing import Optional
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel
from dotenv import load_dotenv
from rdkit import Chem

# Make training/ available on sys.path and import make_diagram
import sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from training.arrowpushing.generation import make_diagram 

load_dotenv()

from openai import OpenAI
client = OpenAI(api_key=os.getenv("VITE_OPENAI_API_KEY") or os.getenv("OPENAI_API_KEY"))

app = FastAPI(title="Chem Chat API", version="1.0")

ALLOWED_ORIGINS = [
    "https://orgomap.org",
    "https://www.orgomap.org",
    # keep local/dev while testing; remove after
    "http://127.0.0.1:5173",
    "http://localhost:5173",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=ALLOWED_ORIGINS,
    allow_credentials=True,
    allow_methods=["POST", "GET", "OPTIONS"],
    allow_headers=["*"],
)

# Static image serving
IMAGES_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "reaction_images"))
os.makedirs(IMAGES_DIR, exist_ok=True)
app.mount("/images", StaticFiles(directory=IMAGES_DIR), name="images")

# --- Models ---
class ChatRequest(BaseModel):
    message: str
    allow_diagrams: bool = True

class ChatResponse(BaseModel):
    mode: str            # "diagram" or "qa"
    text: str            # assistant text answer
    image_url: Optional[str] = None  # when a diagram is generated
    smiles: Optional[str] = None     # reaction smiles used (if any)

# --- Utilities ---
MECH_KEYWORDS = re.compile(
    r"(mechanism|arrow\s*push|curved\s*arrow|show\s*arrows|draw\s*(the\s*)?mechanism|reaction\s*pathway|electron\s*flow)",
    re.IGNORECASE,
)
SMILES_ARROW = re.compile(r">>")

def should_generate_diagram(msg: str) -> bool:
    if MECH_KEYWORDS.search(msg):
        return True
    if SMILES_ARROW.search(msg):
        return True
    if any(k in msg.lower() for k in ["draw", "diagram", "mechanism for", "show the reaction of", "arrow pushing"]):
        return True
    return False

# Minimal local mapper; expand later or hook to LLM if you want
def query_to_smiles(msg: str) -> Optional[str]:
    m = msg.lower().strip()
    def has(*words): return all(w in m for w in words)

    # --- Core presets ---
    # HBr addition to propene (already working)
    if ("propene" in m or "propylene" in m) and ("hbr" in m or "hydrogen bromide" in m or has("h","br")):
        return "C=CC.Br>>CC(Br)C"

    # Nitration of benzene (EAS, HNO3/H2SO4)
    if ("nitration" in m or "nitrate" in m) and ("benzene" in m or "c6h6" in m):
        # agents included as reactants; product = nitrobenzene
        return "c1ccccc1.O=N(=O)O.OS(=O)(=O)O>>O=[N+]([O-])c1ccccc1"

    # Nitration of bromobenzene (EAS) – para/ortho mix; map to p-nitro for demo
    if ("nitration" in m or "nitrate" in m) and "bromobenzene" in m:
        return "Brc1ccccc1.O=N(=O)O.OS(=O)(=O)O>>Brc1ccc([N+](=O)[O-])cc1"

    # Hydration via H2SO4 (ethyl hydrogen sulfate formation from ethene)
    if (("ethene" in m or "ethylene" in m) and ("h2so4" in m or "sulfuric acid" in m)) and \
       any(k in m for k in ["addition", "hydration", "mechanism"]):
        # Captures initial addition to give ethyl hydrogen sulfate
        return "C=C.OS(=O)(=O)O>>CC-OS(=O)(=O)O"  # schematic; RDKit-safe minimal

    # Bromination of alkene (Br2 addition) – test another class
    if any(k in m for k in ["bromination", "add br2", "br2 addition"]) and any(k in m for k in ["alkene", "propene", "ethylene", "cyclohexene"]):
        if "cyclohexene" in m:
            return "C1C=CCCC1.BrBr>>C1C(Br)C(Br)CCC1"
        if "propene" in m or "propylene" in m:
            return "C=CC.BrBr>>C(C)(C)BrBr"  # simplified vic-dibromide sketch
        return "C=C.BrBr>>C(Br)C(Br)"  # generic

    # If user pasted reaction SMILES already, allow
    if ">>" in m:
        return m

    return None



def answer_chem_qa(prompt: str) -> str:
    sysmsg = (
        "You are a concise, precise chemistry tutor."
        " Prefer short, direct answers with necessary equations or structures described in text."
        " If safety or lab procedures are implied, include a brief caution."
    )
    resp = client.chat.completions.create(
        model="gpt-4o-mini",
        temperature=0.2,
        messages=[
            {"role": "system", "content": sysmsg},
            {"role": "user", "content": prompt},
        ],
    )
    return resp.choices[0].message.content.strip()

def safe_filename(seed: str) -> str:
    h = hashlib.sha256(seed.encode("utf-8")).hexdigest()[:16]
    return f"rxn_{h}.png"

@app.post("/api/chemchat", response_model=ChatResponse)
def chemchat(req: ChatRequest):
    msg = req.message.strip()
    if not msg:
        raise HTTPException(400, "Empty message")

    do_diag = req.allow_diagrams and should_generate_diagram(msg)

    smiles: Optional[str] = None
    image_url: Optional[str] = None
    text_answer: str = ""

    try:
        if do_diag:
            if ">>" in msg:
                smiles = msg
            else:
                smiles = query_to_smiles(msg)

            # Validate quickly
            if smiles:
                try:
                    r, p = smiles.split(">>")
                    for part in (r, p):
                        for c in part.split('.'):
                            if c and Chem.MolFromSmiles(c) is None:
                                raise ValueError("Invalid SMILES component")
                except Exception:
                    smiles = None

            if not smiles:
                do_diag = False
                text_answer = (
                    "I couldn't confidently parse a reaction from that. "
                    "Try rephrasing or provide reactants like 'HBr + propene'."
                )

        if do_diag and smiles:
            tmp_svg = os.path.join(IMAGES_DIR, "_tmp.svg")
            out_png = os.path.join(IMAGES_DIR, safe_filename(msg + "::" + smiles))

            # Render SVG (set show_arrows=True when your arrow code is wired)
            make_diagram(smiles, tmp_svg, show_arrows=False)

            # Convert to PNG
            import cairosvg
            cairosvg.svg2png(url=tmp_svg, write_to=out_png)
            if os.path.exists(tmp_svg):
                os.remove(tmp_svg)

            image_url = f"/images/{os.path.basename(out_png)}"
            text_answer = f"Here’s the reaction diagram. Reaction SMILES: `{smiles}`"

        if not text_answer:
            text_answer = answer_chem_qa(msg)

        return ChatResponse(
            mode="diagram" if image_url else "qa",
            text=text_answer,
            image_url=image_url,
            smiles=smiles,
        )

    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(500, f"Error: {e}")

@app.get("/api/health")
def health():
    return {"ok": True}
