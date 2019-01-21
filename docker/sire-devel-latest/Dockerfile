
# This image is used to build, package and test the latest version
# of the devel branch of Sire

FROM siremol/sire-build:latest

WORKDIR $HOME

USER $FN_USER

# Update to the latest version
RUN cd Sire && git checkout devel && git pull
 
# Build Sire (takes a long time)
ADD build_sire.sh $HOME/build_sire.sh
RUN scl enable devtoolset-6 ./build_sire.sh

#Cleaning the conda and pip cache before packaging
ADD clean.py $HOME/clean.py
RUN $HOME/sire.app/bin/conda clean -tipsy && \
    rm -rf $HOME/.cache && \
    $HOME/sire.app/bin/python $HOME/clean.py && \
    rm $HOME/clean.py

ENTRYPOINT ["scl", "enable", "devtoolset-6", "bash"]
