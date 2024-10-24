FROM ubuntu:22.04

WORKDIR /src

RUN apt-get update \
    && apt-get install -y wget \
    && apt install make \
    && apt install python3 -y

# Install BaseSpace CLI
RUN wget https://launch.basespace.illumina.com/CLI/1.6.1/amd64-linux/bs -P /src/tools/bs \
    && chmod a+x /src/tools/bs/bs

# Configure BS CLI
ENV PATH="/src/tools/bs:$PATH"

# Install fastp
RUN wget http://opengene.org/fastp/fastp.0.23.4 -P /src/tools/fastp \
    && mv /src/tools/fastp/fastp.0.23.4 /src/tools/fastp/fastp \
    && chmod a+x /src/tools/fastp/fastp

# Configure fastp
ENV PATH="/src/tools/fastp:$PATH"