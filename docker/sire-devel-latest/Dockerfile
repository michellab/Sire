
# This image is used to build, package and test the latest version
# of the devel branch of Sire

FROM siremol/sire-build:latest

WORKDIR $HOME/Sire

USER $FN_USER

# Update to the latest version
RUN git checkout devel && git pull
 
# Build Sire (takes a long time)
RUN ./compile_sire.sh

ENTRYPOINT ["bash"]
