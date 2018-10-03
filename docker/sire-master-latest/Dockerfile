
# This image is used to build, package and test the latest version
# of the master branch of Sire

FROM chryswoods/sire-build:latest

WORKDIR $HOME

USER $FN_USER

# Update to the latest version
RUN cd Sire && git checkout master && git pull

# Build Sire (takes a long time)
ADD build_sire.sh $HOME/build_sire.sh
RUN scl enable devtoolset-6 ./build_sire.sh

# Now package up into a self-extracting executable
RUN echo "sire_master_latest_linux.run" | $HOME/sire.app/bin/package_sire
RUN chmod a+x sire_master_latest_linux.run

# Now extract this package and run the tests
RUN echo "$HOME/testsire.app" | sire_master_latest_linux.run

RUN $HOME/testsire.app/bin/sire_test 2>&1 | tee sire_tests.log

ENTRYPOINT ["scl", "enable", "devtoolset-6", "bash"]
