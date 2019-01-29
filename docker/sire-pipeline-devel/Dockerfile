
# This image is used to build, package and test the latest version
# of the devel branch of Sire

FROM siremol/sire-devel:latest

WORKDIR $HOME

USER $FN_USER

# Update to the latest version
RUN cd Sire && git checkout devel && git pull
 
# Build Sire (this should take less time as it is a small update)
RUN scl enable devtoolset-6 ./build_sire.sh

ENTRYPOINT ["scl", "enable", "devtoolset-6", "bash"]
