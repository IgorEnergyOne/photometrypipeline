#!/bin/bash

# Photometrypipeline installation script
# Based on https://photometrypipeline.readthedocs.io/en/latest/install.html

# Install package dependencies
sudo apt-get update
sudo apt-get install -y \
       build-essential \
       libssl-dev \
       libffi-dev \
       git \
       wget \
       imagemagick \
       curl \
       libplplot-dev \
       libshp-dev \
       libcurl4-gnutls-dev \
       liblapack3 liblapack-dev liblapacke liblapacke-dev \
       libfftw3-3 libfftw3-dev libfftw3-single3 \
       libatlas-base-dev \
       sextractor \
       scamp

### Install Python dependencies
yes | pip install -r requirements.txt


# check if system is Linux Mint distribution
if grep -q "NAME=\"Linux Mint\"" /etc/os-release; then
  echo "Linux Mint detected. Adding photometrypipeline to PATH"
  # photometry pipeline setup
  # add photometrypipeline to path
  echo 'export PATH="$PATH:'"$(pwd)"'"' >> ~/.bashrc
  echo "export PHOTPIPEDIR='$(pwd)'" >> ~/.bashrc
else
  echo "Please add photometrypipeline to PATH manually:"
  echo "#photometry pipeline setup"
  echo "export PHOTPIPEDIR=<path>/photometrypipeline"
  echo "export PATH=$PATH:~<path>/photometrypipeline/"
fi