# Image reproductible de msannot : versions exactes (requirements.txt), utilisateur non root.
#   docker build -t msannot .
#   docker run --rm -v "$PWD/results/docker:/app/results" msannot           # démonstration
#   docker run --rm -p 8501:8501 --entrypoint streamlit msannot \
#       run app/streamlit_app.py --server.address 0.0.0.0                   # dashboard
FROM python:3.12-slim

# Bibliothèques X11 requises par le rendu des molécules de RDKit (Cairo).
RUN apt-get update \
    && apt-get install -y --no-install-recommends libxrender1 libxext6 \
    && rm -rf /var/lib/apt/lists/*

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    MPLBACKEND=Agg \
    MPLCONFIGDIR=/tmp/matplotlib

WORKDIR /app
COPY requirements.txt ./
RUN pip install -r requirements.txt

COPY pyproject.toml README.md LICENSE ./
COPY src ./src
RUN pip install --no-deps .

COPY config ./config
COPY data/demo ./data/demo
COPY app ./app

RUN useradd --create-home --uid 1000 msannot \
    && mkdir -p /app/results \
    && chown msannot /app/results
USER msannot

EXPOSE 8501
ENTRYPOINT ["msannot"]
CMD ["demo", "--output", "/app/results/demo"]
