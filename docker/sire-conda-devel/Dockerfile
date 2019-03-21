# This image is used to build the Sire Conda package.

FROM siremol/sire-devel:latest

WORKDIR $HOME

USER $FN_USER

# Import arguments and set environment variables.
ARG par_url
ARG github_token
ARG github_email
ENV PAR_URL=$par_url
ENV GITHUB_TOKEN=$github_token
ENV GITHUB_EMAIL=$github_email

# Create the compressed Conda package.
ADD package_conda.sh .
RUN $HOME/package_conda.sh

# Install pycurl
RUN $HOME/sire.app/bin/conda install -y pycurl

# Deploy the Conda package.
ADD deploy.py .
RUN $HOME/sire.app/bin/python deploy.py $HOME/sire_conda_*.tar.bz2

# Udpate the Conda Forge feedstock.
ADD update_feedstock.sh
RUN $HOME/update_feedstock.sh $GITHUB_TOKEN $GITHUB_EMAIL

ENTRYPOINT ["bash"]
