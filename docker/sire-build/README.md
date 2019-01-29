How to run the docker Docker file contained here:
-----------------------
`docker build -t siremol/sire-build:latest .`


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
`docker run -it siremol/sire-build:latest --netowrk=host`

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
