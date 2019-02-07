
# This image is used to build, package and test the latest version
# of the devel branch of Sire

FROM siremol/sire-build:latest

WORKDIR $HOME/Sire

USER $FN_USER

# Update to the latest version
RUN git checkout devel && git pull
 
# Build Sire (takes a long time)
RUN ./compile_sire.sh

#Cleaning the conda and pip cache before packaging
ADD clean.py $HOME/clean.py
RUN $HOME/sire.app/bin/conda clean -tipsy && \
    rm -rf $HOME/.cache && \
    $HOME/sire.app/bin/python $HOME/clean.py && \
    rm $HOME/clean.py

ENTRYPOINT ["bash"]
