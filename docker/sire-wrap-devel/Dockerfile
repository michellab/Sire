# This image is used to create a container for building Sire Python wrappers.

FROM siremol/sire-devel:latest

# Configure environment
ENV SHELL=/bin/bash \
    FN_USER=sireuser \
    FN_UID=1000 \
    FN_GID=100 \
    SIRE_SILENT_PHONEHOME=1 \
    SIRE_SILENT_DONT_PHONEHOME=1 \
    HOME=/home/$FN_USER \
    PATH=/home/$FN_USER/sire.app/bin:$PATH \
    LD_LIBRARY_PATH=$HOME/$FN_USER/sire.app/lib:$HOME/$FN_USER/sire.app/lib64:$LD_LIBRARY_PATH

USER $FN_USER

# Install pyplusplus and pygccxml.
RUN pip install pyplusplus fuzzywuzzy && \
    conda install -y -c zeus1942 pygccxml && \
    conda install -y -c conda-forge clang clangdev llvmdev gcc_linux-64 gxx_linux-64

# Download and install CastXML.
RUN cd $HOME && \
    git clone https://github.com/CastXML/CastXML.git && \
    cd CastXML && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_INSTALL_PREFIX=$HOME/sire.app .. && \
    make -j 4 && \
    make install

# Clean up so that the Docker image size is minimised.
RUN conda clean -tipy

# Set the entry directory.
WORKDIR $HOME/Sire/wrapper

ENTRYPOINT ["bash"]
