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

# Clone HiGHS repository
RUN git clone https://github.com/ERGO-Code/HiGHS.git

# Build HiGHS
WORKDIR /HiGHS/src
RUN mkdir build && \
    cd build && \
    cmake .. && \
    cmake --build . && \
    ctest

WORKDIR /pypsa-earth-trace

ENV PATH="/opt/conda/envs/base/bin:/HiGHS/src/build/bin:$PATH"

ARG MAMBA_DOCKERFILE_ACTIVATE=1  # (otherwise python will not be found)

RUN python -c 'import uuid; print(uuid.uuid4())' > /tmp/my_uuid

COPY . /pypsa-earth-trace