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

# Install fastqc (will include dependencies like java)
RUN apt install fastqc=0.11.9+dfsg-5 -y

# Install STAR aligner
RUN wget https://github.com/alexdobin/STAR/archive/2.7.11b.tar.gz -P /src/tools/ \
    && tar xzvf /src/tools/2.7.11b.tar.gz -C ./tools \
    && rm /src/tools/2.7.11b.tar.gz

# Configure STAR alginer
ENV PATH="/src/tools/STAR-2.7.11b/bin/Linux_x86_64_static:$PATH"

# Install pip
RUN apt install python3-pip -y

# Install pydeseq2
RUN pip3 install pydeseq2==0.4.12

# Install scanpy
RUN pip3 install scanpy==1.10.3

# Install pyComBat via InMoose
RUN pip3 install --upgrade pip
RUN pip3 install inmoose --break-system-packages