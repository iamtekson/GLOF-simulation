#!/bin/bash

#Making sure that Python 3 is used
sudo apt install python3-pip -y
echo "alias python=python3" >> ~/.bash_aliases
source ~/.bash_aliases

#Selecting suitable repository
sudo add-apt-repository ppa:ubuntugis/ubuntugis-unstable -y
sudo apt update

#Installing necessary additional packages
sudo apt install libgeos-dev -y
sudo apt install libproj-dev proj-data proj-bin -y
sudo apt install libgdal-dev python3-gdal gdal-bin -y

#Installing GRASS GIS (dev package)
sudo apt update
sudo apt install grass-dev grass-doc grass-gui -y

#Installing pillow
sudo apt install python3-pip
python3 -m pip install --upgrade pillow

#Installing R statistical software
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
sudo apt update
sudo apt install r-base r-base-core -y

#Installing necessary additional R packages
echo 'install.packages(c("stats","foreign","sp","rgeos","rgdal","raster","maptools","ROCR","fmsb", "Rcpp"), lib="/usr/lib/R/library/", repos = "http://cran.case.edu" )' > install.packages.R
sudo R CMD BATCH install.packages.R

#Terminating sudo mode
sudo -k
#clear
