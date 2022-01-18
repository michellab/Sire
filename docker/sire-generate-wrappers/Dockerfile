FROM alpine:latest

RUN apt-get update && apt-get -y upgrade \
  && apt-get install -y --no-install-recommends \
    git \
    wget \
    g++ \
    gcc \
    nano \
    ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh \
    && echo "Running $(conda --version)" && \
    conda init bash && \
    . /root/.bashrc && \
    conda install make clang clangdev llvmdev cmake && \
    conda clean -a -f -y && \
    pip install pyplusplus pygccxml fuzzywuzzy python-Levenshtein && \
    rm -fr ~/.cache/pip /tmp*

ENV LD_LIBRARY_PATH="/root/miniconda3/lib:/root/miniconda3/lib64:${LD_LIBRARY_PATH}"
WORKDIR /root

RUN git clone https://github.com/CastXML/CastXML && \
    cd CastXML && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_INSTALL_PREFIX=$HOME/miniconda3/ .. && \
    make -j 4 && \
    make -j 4 install && \
    cd $HOME && \
    rm -rf CastXML

COPY includes.tar.bz2 /root/miniconda3/
COPY prepare.sh /root/

