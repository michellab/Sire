
# This image is used to build, package and test the latest version
# of the devel branch of Sire

FROM siremol/sire-package-devel:latest

WORKDIR $HOME

USER $FN_USER

ARG par_url
ENV PAR_URL=$par_url

ADD deploy.py .
RUN $HOME/sire.app/bin/python deploy.py

ENTRYPOINT ["scl", "enable", "devtoolset-6", "bash"]
