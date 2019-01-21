
# This image is used to build, package and test the latest version
# of the devel branch of Sire

FROM siremol/sire-devel:latest

WORKDIR $HOME

USER $FN_USER
#Cleaning the conda and pip cache before packaging
RUN $HOME/sire.app/bin/conda clean -tipsy
RUN rm -rf $HOME/.cache

# Now package up into a self-extracting executable
RUN echo "sire_devel_latest_linux.run" | $HOME/sire.app/bin/package_sire
RUN chmod a+x sire_devel_latest_linux.run

ENTRYPOINT ["scl", "enable", "devtoolset-6", "bash"]
