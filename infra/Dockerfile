# syntax = docker/dockerfile:1
FROM mambaorg/micromamba:bookworm-slim AS build-stage

USER root

ENV PIP_NO_CACHE_DIR=1
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

RUN apt-get update && apt-get install -y --no-install-recommends \
        git \
        cmake \
        g++ \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN --mount=type=bind,source=envs/environment.yaml,target=/tmp/env.yaml \
    micromamba create -y --always-copy -p /env -f /tmp/env.yaml \
    && micromamba clean --all --force-pkgs-dirs --yes \
    && rm -rf /env/conda-meta \
    && rm -rf /env/include \
    && find -name '*.a' -delete \
    && find -name '__pycache__' -type d -exec rm -rf '{}' '+'

FROM gcr.io/distroless/base-debian12 AS run-stage

WORKDIR /pypsa-earth-trace

COPY --from=build-stage /env /env
ENV PATH="/env/bin:${PATH}"

COPY . .
