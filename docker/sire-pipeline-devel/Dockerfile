
# This image is used to build, package and test the latest version
# of the devel branch of Sire

FROM siremol/sire-devel:latest

WORKDIR $HOME/Sire

USER $FN_USER

# Update to the latest version
RUN git checkout devel && git pull
 
# Build Sire (this should take less time as it is a small update)
RUN ./compile_sire.sh

ENTRYPOINT ["bash"]
