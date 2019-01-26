
# This image is used to build, package and test the latest version
# of the devel branch of Sire

FROM siremol/sire-devel:latest

WORKDIR $HOME

USER $FN_USER

RUN $HOME/sire.app/bin/sire_test 2>&1 | tee sire_tests.log

ENTRYPOINT ["scl", "enable", "devtoolset-6", "bash"]
