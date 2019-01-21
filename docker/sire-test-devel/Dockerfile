
# This image is used to build, package and test the latest version
# of the devel branch of Sire

FROM siremol/sire-package:latest

WORKDIR $HOME

USER $FN_USER

# Now extract this package and run the tests
RUN echo "$HOME/testsire.app" | ./sire_devel_latest_linux.run

RUN $HOME/testsire.app/bin/sire_test 2>&1 | tee sire_tests.log

ENTRYPOINT ["scl", "enable", "devtoolset-6", "bash"]
