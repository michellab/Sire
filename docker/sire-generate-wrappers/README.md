# Creating the Sire wrappers docker container

Aim is to create a Docker container that has a working
Py++ and nothing else (i.e. not a full Sire install).

Discovered that we need to also include the header files
for the dependencies of Sire. Can get these by 
tarring them up from a separate successful Sire install.

Assuming you have already installed Sire into $HOME/sire.app,
you can create the necessary includes.tar.bz2 by running

```
cd $HOME/sire.app
tar -jcvf /path/to/sire_wrappers_docker/includes.tar.bz2 \
      include/qt include/gsl include/openmm \
      include/OpenMM* include/boost include/tbb \
      include/oneapi
cd -
```

This creates a tarball of all of the needed headers, which
are unpacked within the docker container into the 
miniconda3/include directory.

