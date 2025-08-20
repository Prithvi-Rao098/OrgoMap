# Hugging Face Spaces loves 7860; micromamba + rdkit is lightweight and current
FROM mambaorg/micromamba:1.5.10

ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN micromamba install -y -n base -c conda-forge \
      python=3.10 rdkit=2024.09.* \
    && micromamba clean -a -y

# libs for cairosvg PNG rendering
USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
      libcairo2 libpango-1.0-0 fonts-dejavu-core \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY server ./server
COPY training ./training

RUN pip install --no-cache-dir fastapi uvicorn pydantic python-dotenv openai cairosvg

ENV PORT=7860
EXPOSE 7860
CMD uvicorn server.main:app --host 0.0.0.0 --port ${PORT}
