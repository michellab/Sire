# This image is used to build and deploy the Sire Conda package.

FROM siremol/sire-devel:latest

WORKDIR $HOME

USER $FN_USER

# Import arguments and set environment variables.
ARG par_url
ENV PAR_URL=$par_url

# Create the compressed Conda package.
ADD create_package_file.sh .
RUN $HOME/create_package_file.sh

# Install PycURL, conda-build, and anaconda-client
RUN $HOME/sire.app/bin/conda install -y pycurl conda-build anaconda-client

# Deploy the Conda package files.
ADD deploy.py .
RUN $HOME/sire.app/bin/python deploy.py $HOME/sire_conda_latest_linux.tar.bz2

# Update the Conda package recipe.
ADD update_recipe.sh .
RUN $HOME/update_recipe.sh

ENTRYPOINT ["bash"]
