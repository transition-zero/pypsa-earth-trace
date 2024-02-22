FROM mambaorg/micromamba

WORKDIR /pypsa-earth-trace

COPY --chown=$MAMBA_USER:$MAMBA_USER envs/environment.yaml /tmp/env.yaml

USER root

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y git cmake g++ && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set the GRB_LICENSE_FILE environment variable
ENV GRB_LICENSE_FILE=/mnt/gcs/gurobi.lic

RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes



WORKDIR /pypsa-earth-trace

ARG MAMBA_DOCKERFILE_ACTIVATE=1  # (otherwise python will not be found)

RUN python -c 'import uuid; print(uuid.uuid4())' > /tmp/my_uuid

COPY . /pypsa-earth-trace