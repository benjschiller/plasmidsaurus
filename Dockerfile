FROM python:3.12.7-bullseye AS base

# Versions
ARG SPADES_VERSION=4.0.0
ARG PROKKA_VERSION=v1.14.5

ENV SPADES_VERSION=${SPADES_VERSION}
ENV PROKKA_VERSION=${PROKKA_VERSION}

# Install OS-level dependencies
RUN apt-get update && apt-get install -y \
    wget \
    # these are for prokka \
    libdatetime-perl libxml-simple-perl libdigest-md5-perl git default-jre bioperl \
    && rm -rf /var/lib/apt/lists/*

# Install BioPerl
RUN cpan Bio::Perl

# Make /opt directory for eveything
RUN mkdir -p /opt
WORKDIR /opt

# Install spades in a separate layer (debug-friendly)
FROM base AS spades
RUN wget -q https://github.com/ablab/spades/releases/download/v${SPADES_VERSION}/SPAdes-${SPADES_VERSION}-Linux.tar.gz && \
    tar xzf SPAdes-${SPADES_VERSION}-Linux.tar.gz -C /opt && \
    rm -f /opt/SPAdes-${SPADES_VERSION}-Linux.tar.gz


# Install prokka in a spearate layer
FROM base AS prokka
RUN wget -q https://github.com/tseemann/prokka/archive/refs/tags/${PROKKA_VERSION}.tar.gz && \
    tar xzf ${PROKKA_VERSION}.tar.gz && \
    # n.b. inconsistent use of v in naming, add it back
    mv prokka-${PROKKA_VERSION#v} prokka-${PROKKA_VERSION} && \
    rm -f ${PROKKA_VERSION}.tar.gz
# expose included binaries in PATH (smdh)
ENV PATH="/opt/prokka-${PROKKA_VERSION}/binaries/linux:$PATH"
RUN /opt/prokka-${PROKKA_VERSION}/bin/prokka --setupdb

# Install venv separately
FROM base AS venv
# Create a virtual environment
RUN python3 -m venv venv
ADD requirements.txt /opt/requirements.txt

# Install the Python dependencies
RUN /opt/venv/bin/pip install --no-cache-dir -r requirements.txt

# bring it all together
FROM base
COPY --from=spades /opt/SPAdes-${SPADES_VERSION}-Linux /opt/SPAdes-${SPADES_VERSION}-Linux
COPY --from=prokka /opt/prokka-${PROKKA_VERSION} /opt/prokka-${PROKKA_VERSION}
COPY --from=venv /opt/venv /opt/venv

# Set the environment variable to use the virtual environment
ENV VIRTUAL_ENV=/opt/venv
ENV PATH="/opt/venv/bin:/opt/prokka-${PROKKA_VERSION}/bin:/opt/prokka-${PROKKA_VERSION}/binaries/common:/opt/prokka-${PROKKA_VERSION}/binaries/linux:/opt/SPAdes-${SPADES_VERSION}-Linux/bin:$PATH"
