
# This image is used to build, package and test the latest version
# of the master branch of Sire

FROM chryswoods/sire-build:latest

WORKDIR $HOME

USER $FN_USER

# Update to the latest version
RUN cd Sire && git checkout master && git pull

ENTRYPOINT ["bash"]
