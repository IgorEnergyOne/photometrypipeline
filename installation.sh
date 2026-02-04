#!/bin/bash

# Photometrypipeline installation script

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# 1. Check for Anaconda/Miniconda
if command_exists conda; then
    echo "Conda detected."
else
    echo "Conda not found. Installing Miniconda..."
    mkdir -p ~/miniconda3
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
    bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
    rm -rf ~/miniconda3/miniconda.sh

    # Initialize conda for bash
    ~/miniconda3/bin/conda init bash
    source ~/.bashrc
fi

# 2. Install package dependencies via apt-get
echo "Installing system dependencies..."
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

# 3. Check Python in base environment
echo "Checking Python in Conda base environment..."
# Check if python is running from conda
current_python=$(which python)
if [[ $current_python == *"conda"* ]]; then
    echo "Python is already installed in the current Conda environment."
else
    echo "Ensuring Python is installed in base..."
    # Attempt to activate base just in case, though script context might differ
    # Installing python 3.10 as a safe default for modern scientific stacks
    conda install -n base python=3.10 -y
fi

# 4. Install Python dependencies
echo "Installing Python requirements..."
pip install -r requirements.txt

# 5. Configure PATH based on OS
OS_NAME=$(grep ^NAME= /etc/os-release | cut -d= -f2 | tr -d '"')
IS_WSL=$(grep -i microsoft /proc/version)
CURRENT_DIR=$(pwd)

echo "Configuring environment variables..."

# Detect if WSL
if [[ -n "$IS_WSL" ]]; then
    echo "WSL detected."
    # WSL often shares .bashrc standard behavior, but specific adjustments can be made here if needed.
    # Appending to .bashrc works standardly for WSL bash.
fi

if [[ "$OS_NAME" == "Linux Mint" ]] || [[ "$OS_NAME" == "Ubuntu" ]] || [[ -n "$IS_WSL" ]]; then
    echo "Compatible distribution detected ($OS_NAME). Updating .bashrc..."

    # Avoid duplicate entries
    if ! grep -q "export PHOTPIPEDIR='$CURRENT_DIR'" ~/.bashrc; then
        echo "" >> ~/.bashrc
        echo "# photometry pipeline setup" >> ~/.bashrc
        echo "export PHOTPIPEDIR='$CURRENT_DIR'" >> ~/.bashrc
        echo 'export PATH="$PATH:$PHOTPIPEDIR"' >> ~/.bashrc
        echo "Path updated. Please restart your terminal or run 'source ~/.bashrc'."
    else
        echo "Photometry pipeline already in PATH."
    fi
else
    echo "Distribution not explicitly matched for auto-config ($OS_NAME)."
    echo "Please add the following manually to your shell configuration:"
    echo "export PHOTPIPEDIR='$CURRENT_DIR'"
    echo 'export PATH="$PATH:$PHOTPIPEDIR"'
fi
