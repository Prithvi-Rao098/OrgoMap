# OrgoMap

OrgoMap is an **AI-powered organic chemistry tool** that generates **reaction diagrams** (with optional arrow-pushing mechanisms) directly from **reaction SMILES**.

This repo contains the **backend service** for the beta release, built with **FastAPI**, **RDKit**, and deployed on **Hugging Face Spaces (Docker runtime).**

---

## 🚀 Features
- Parse **reaction SMILES → RDKit objects**  
- Generate **clean reaction diagrams (SVG)** without atom indices  
- (Beta) Support for **arrow-pushing electron flow**  
- Batch-generate diagrams from **USPTO parquet datasets**  
- Expose a **REST API** for frontend integration  

---

## ⚡ API Endpoints
Once deployed, your backend will be hosted at:

https://<username>-orgomap-backend.hf.space

Available routes:
- `GET /` → Health check (`{"Hello": "World!"}`)  
- `POST /diagram` → Generate a single diagram from SMILES  
- `POST /batch` → Batch-generate diagrams  

---

## 🧑‍💻 Local Setup

### 1. Clone the repo
```bash
git clone https://huggingface.co/spaces/PrithviRao/OrgoMap-backend
cd OrgoMap-backend
2. Install dependencies
conda create -n orgomap python=3.9 -y
conda activate orgomap
pip install -r requirements.txt
3. Run server
uvicorn app:app --reload --host 0.0.0.0 --port 7860
Visit → http://localhost:7860

🛠 Deployment (Hugging Face Spaces)
Push this repo to your Hugging Face Space

In Settings → Variables & secrets, add:
OPENAI_API_KEY = your_real_key
Ensure Runtime = Docker

HF will build automatically → backend live at:

https://<username>-orgomap-backend.hf.space
🧪 Example Usage
from training.arrowpushing import generation

# Single diagram
generation.make_diagram(
    "CC(C)(C)Br.[OH-]>>C=C(C)C.Br.[H2O]",
    out_svg="elim.svg",
    show_arrows=True
)

# Batch diagrams
generation.batch_generate_diagrams(
    "uspto.parquet",
    "./out_diagrams"
)
Output → elim.svg

🔬 Chemistry Example
E2 Elimination (tert-butyl bromide + OH⁻ → isobutylene):
CC(C)(C)Br.[OH-]>>C=C(C)C.Br.[H2O]
⚠️ Notes
Hugging Face blocks large binaries (PNG, parquet)

Only lightweight training code + examples included

Arrow-pushing support is still experimental

📜 License
MIT License – free to use, modify, and distribute
