
# This image is used to build, package and test the latest version
# of the devel branch of Sire

FROM siremol/sire-build:latest

WORKDIR $HOME/Sire

USER $FN_USER

# Update to the latest version
RUN git checkout devel && git pull
 
# Build Sire (takes a long time)
RUN ./compile_sire.sh

# Clean up so that the Docker image size is minimised
RUN $HOME/sire.app/bin/conda clean -tipy

ENTRYPOINT ["bash"]
