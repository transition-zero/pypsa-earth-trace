FROM mambaorg/micromamba

WORKDIR /pypsa-earth-trace

COPY --chown=$MAMBA_USER:$MAMBA_USER envs/environment.yaml /tmp/env.yaml

USER root

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y git cmake g++ && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1  # (otherwise python will not be found)

COPY . .
