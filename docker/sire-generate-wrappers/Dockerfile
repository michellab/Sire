# This image is used to generate and commit Sire Python wrappers.

FROM sire/sire-wrap-devel:latest

WORKDIR $HOME

ARG github_token
ARG github_email
ENV GITHUB_TOKEN=$github_token
ENV GITHUB_EMAIL=$github_email

# Update the Sire source code.
RUN git pull

# Generate the wrappers.
RUN $HOME/sire.app/bin/python AutoGenerate/scanheaders.py $HOME/Sire/corelib/src/libs . && \
    $HOME/sire.app/bin/python create_all_wrappers.py

# Commit and push any new/updated Python wrappers.
RUN git config user.name "BioSimSpaceBot" && \
    git config user.email "$GITHUB_EMAIL" && \
    git commit -a -m "Updating Python wrappers." && \
    git push --repo https://biosimspacebot:$GITHUB_TOKEN@github.com/michellab/Sire.git > /dev/null 2>&1

ENTRYPOINT ["bash"]
