FROM continuumio/miniconda3:latest

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

WORKDIR /root

RUN conda install make clang clangdev llvmdev cmake && \
    conda clean -a -f -y && \
    pip install pyplusplus pygccxml fuzzywuzzy python-Levenshtein && \
    rm -fr ~/.cache/pip /tmp*

RUN git clone https://github.com/CastXML/CastXML && \
    cd CastXML && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_INSTALL_PREFIX=/opt/conda/ .. && \
    make -j 4 && \
    make -j 4 install && \
    cd $HOME && \
    rm -rf CastXML

COPY includes.tar.bz2 /opt/conda
COPY generate_wrappers /usr/bin
RUN chmod a+x /usr/bin/generate_wrappers
COPY bashrc /root/.bashrc
COPY push_wrappers /usr/bin
RUN chmod a+x /usr/bin/push_wrappers

RUN mkdir /tmp
