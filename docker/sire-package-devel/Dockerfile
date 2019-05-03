
# This image is used to build, package and test the latest version
# of the devel branch of Sire

FROM siremol/sire-devel:latest

WORKDIR $HOME

USER $FN_USER

#Cleaning the conda and pip cache before packaging
RUN $HOME/sire.app/bin/conda uninstall --yes cmake make gcc_linux-64 gxx_linux-64
RUN $HOME/sire.app/bin/conda clean -tipy
RUN rm -rf $HOME/.cache

# Now package up into a self-extracting executable
ENV SIRE_RUN_FILE sire_devel_latest_linux.run
RUN $HOME/sire.app/bin/package_sire > package_sire.log 2> package_sire.err
RUN chmod a+x $SIRE_RUN_FILE

ENTRYPOINT ["bash"]
