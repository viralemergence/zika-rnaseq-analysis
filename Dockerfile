FROM ubuntu:22.04

WORKDIR /src

RUN apt-get update \
    && apt-get install -y wget

# Install BaseSpace CLI
RUN wget https://launch.basespace.illumina.com/CLI/latest/amd64-linux/bs -P /src/tools/bs \
    && chmod u+x /src/tools/bs/bs

# Configure BS CLI
ENV PATH="/src/tools/bs:$PATH"