
# This image is used to build, package and test the latest version
# of the devel branch of Sire

FROM siremol/sire-package:latest

WORKDIR $HOME

USER $FN_USER

# Set installation directory.
ENV SIRE_INSTALL_DIR $HOME/test_sire.app

RUN ./sire_devel_latest_linux.run

RUN $SIRE_INSTALL_DIR/bin/sire_test 2>&1 | tee sire_tests.log

ENTRYPOINT ["scl", "enable", "devtoolset-6", "bash"]
