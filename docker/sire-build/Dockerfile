# We build Sire in Scientific Linux 6. This has an old glibc that
# means that binaries are compatible with a wide array of Linux distros

# This image is used to create the base container for all Sire builds

FROM sl:6
ENV container docker

RUN yum -y install git wget

# Configure environment
ENV SHELL=/bin/bash \
    FN_USER=sireuser \
    FN_UID=1000 \
    FN_GID=100

ENV HOME=/home/$FN_USER

ADD fix-permissions /usr/bin/fix-permissions

# Create fn user with UID=1000 and in the 'users' group
# and make sure these dirs are writable by the `users` group.
RUN useradd -m -s /bin/bash -N -u $FN_UID $FN_USER && \
    fix-permissions $HOME

WORKDIR $HOME

RUN /usr/bin/fix-permissions $HOME

USER $FN_USER

ENV PATH=$HOME/local/bin:$PATH

# Seed the repo with the latest devel version of Sire.
# This will stop us downloading too much for each build...
RUN git clone https://github.com/michellab/Sire

ENTRYPOINT ["bash"]
