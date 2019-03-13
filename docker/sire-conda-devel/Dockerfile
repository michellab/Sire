# This image is used to build the Sire Conda package.

FROM siremol/sire-devel:latest

WORKDIR $HOME

USER $FN_USER

# Import arguments and set environment variables.
#ARG par_url
#ENV PAR_URL=$par_url

# Add the packaging script.
ADD package_conda.sh .

# Create the compressed Conda package.
RUN $HOME/package_conda.sh

# Install pycurl
#RUN $HOME/sire.app/bin/conda install -y pycurl

# Deploy the Conda package.
#ADD deploy.py .
#RUN $HOME/sire.app/bin/python deploy.py $HOME/sire_conda_*.tar.bz2

ENTRYPOINT ["bash"]
