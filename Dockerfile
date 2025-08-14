# Dockerfile
FROM rdkit/minimal:latest

WORKDIR /app
# copy backend and RDKit code
COPY server ./server
COPY training ./training

# python deps
RUN pip install --no-cache-dir fastapi uvicorn pydantic python-dotenv openai cairosvg

# images dir (ephemeral on Render, fine for demo)
RUN mkdir -p /app/reaction_images
ENV PORT=10000

CMD uvicorn server.main:app --host 0.0.0.0 --port $PORT
