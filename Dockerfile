FROM python:3.12.7-bullseye AS base

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Versions
ARG MINIFORGE_VERSION="24.7.1-2"
ARG FILTLONG_VERSION=v0.2.1
ARG FLYE_VERSION=2.9.5
ARG PROKKA_VERSION=v1.14.5

ENV MINIFORGE_VERSION=${MINIFORGE_VERSION}
ENV FILTLONG_VERSION=${FILTLONG_VERSION}
ENV PROKKA_VERSION=${PROKKA_VERSION}

# Install OS-level dependencies
RUN apt-get update && apt-get install -y \
    wget git \
    && rm -rf /var/lib/apt/lists/*

# Make /opt directory for eveything
RUN mkdir -p /opt
WORKDIR /opt

FROM base AS mamba
# Install Miniforge3 (mamba and conda)
RUN curl -L -O "https://github.com/conda-forge/miniforge/releases/download/${MINIFORGE_VERSION}/Miniforge3-${MINIFORGE_VERSION}-$(uname)-$(uname -m).sh" && \
    bash Miniforge3-${MINIFORGE_VERSION}-$(uname)-$(uname -m).sh -b -p /opt/miniforge3
COPY environment.yaml .
RUN /opt/miniforge3/bin/conda env create -f environment.yaml

# Install Filtlong (use separate layers for easy debugging and faster builds)
FROM base AS filtlong
RUN wget -q https://github.com/rrwick/Filtlong/archive/refs/tags/${FILTLONG_VERSION}.tar.gz && \
    tar xfz ${FILTLONG_VERSION}.tar.gz && \
    mv Filtlong-${FILTLONG_VERSION#v} Filtlong-${FILTLONG_VERSION} && \
    cd Filtlong-${FILTLONG_VERSION} && \
    make && \
    bin/filtlong -h

# bring it all together
FROM mamba
COPY --from=filtlong /opt/Filtlong-${FILTLONG_VERSION} /opt/Filtlong-${FILTLONG_VERSION}

# Set the default conda env
ENV CONDA_DEFAULT_ENV $conda_env

ENV PATH="/opt/Filtlong-${FILTLONG_VERSION}/bin:$PATH"

COPY entrypoint.sh .
ADD . /opt/app
WORKDIR /opt/app
ENTRYPOINT [ "/opt/entrypoint.sh" ]