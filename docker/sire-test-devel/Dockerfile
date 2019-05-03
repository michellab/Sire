
# This image is used to build, package and test the latest version
# of the devel branch of Sire

FROM siremol/sire-devel:latest

WORKDIR $HOME

USER $FN_USER

RUN $HOME/sire.app/bin/sire_test

ENTRYPOINT ["bash"]
