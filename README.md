# contour2DSO

Creation of a DICOM Segmentation Object (DSO) from an AIM file storing the contours of a DICOM image series. Integrated into Stanford University's ePAD Imaging Platform as a plugin. View full documentation on the [ePAD website](https://epad.stanford.edu/contour-dso).

Available on Docker Hub [here](https://hub.docker.com/repository/docker/rvignav/contour2dso).

Table of Contents
=================

<!--ts-->
   * [Usage](#usage)
   * [Errors](#errors)
<!--te-->

Usage
=====

To test the contour to DSO algorithm, run the following commands:

    git clone https://github.com/rvignav/contour2DSO.git
    cd contour2DSO
    docker build -t run .
    docker run -it -v "local/aims:/home/series/files" -v "local/patient/series:/home/series/PatientSeries" run

A possible command satisfying the bind and argument requirements is:

    docker run -it -v "$(pwd)/aims:/home/series/files" -v "$(pwd)/SamplePatient:/home/series/PatientSeries" run

Note that the AIM file must be stored in the `aims` folder.

The generated DSO is now stored in the `/output` folder of the Docker container and can be accessed by ePAD.

Errors
=====

To fix the error `Cannot connect to the Docker daemon at unix:///var/run/docker.sock. Is the docker daemon running?`, run:

Windows:

    systemctl start docker

MacOS:

    brew cask install docker virtualbox
    brew install docker-machine
    docker-machine create --driver virtualbox default
    docker-machine restart
    eval "$(docker-machine env default)"

To fix the error `docker: Error response from daemon: Conflict. The container name <container-name> is already in use`, run:

    docker ps -q -a | xargs docker rm
