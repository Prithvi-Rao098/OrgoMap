# RDKit base with Python
FROM mcs07/rdkit:latest

# (Optional) Cairo runtime & fonts for cairosvg PNG rendering
RUN apt-get update && apt-get install -y libcairo2 fonts-dejavu-core && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY server ./server
COPY training ./training

# Python deps your API needs
RUN pip install --no-cache-dir fastapi uvicorn pydantic python-dotenv openai cairosvg

# Render/Fly set PORT; default to 10000 if missing
ENV PORT=10000
CMD uvicorn server.main:app --host 0.0.0.0 --port ${PORT}