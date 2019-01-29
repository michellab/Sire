How to run the docker Docker file contained here:
-----------------------
`docker build -t siremol/sire-package:latest .`

For packaging choose whichever sire build you like and adjust the Dockerfile accordingly.    
In order to get the Sire run binary out of the container please see below how to copy data out of a container. 


Docker: a quick intro
---------------------

In order to run the docker container make sure you have docker installed.
Type `docker` in your command line.

### Help:
`docker --help` will provide you with useful information 

### Check available images:
`docker images` gives you a list of available docker images

### Run interactively:
If you want to interactively run a docker file to check output:   
`docker run -it siremol/sire-package:latest --netowrk=host`

### Copy data out of a docker container:   
1. Check your docker contaienr is running with:   
`docker ps` a list of currently running docker processes   
If it isn't running you can start it with the example above interactively. 
2. Find your process ID with `docker ps`:   
`CONTAINER ID      c32263c3eb25`   
3. Copy data using your process ID:
`docker cp pocessID:/full/path/to/file /path/where/you/want/it`
An example would be:   
`docker cp c32263c3eb25:/home/sireuser/compile_sire.log .`
