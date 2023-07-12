#!/bin/bash

# conda env
conda env create -f environment.yml

# git dependencies seqrequester/build/bin/seqrequester
if ! command -v seqrequester &> /dev/null; then
    echo "seqrequester not found in PATH. Installing..."
    git clone https://github.com/marbl/seqrequester.git
    cd seqrequester/src
    make
    cd ../../
    echo "Installation complete!"
else
    echo "Unable to install read simulator"
fi

https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.5/sratoolkit.3.0.5-ubuntu64.tar.gz

# git dependencies seqrequester/build/bin/seqrequester
if ! command -v prefetch &> /dev/null; then
    echo "prefetch not found in PATH. Installing SRA tools..."
    wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.5/sratoolkit.3.0.5-ubuntu64.tar.gz
    cd seqrequester/src
    make
    cd ../../
    echo "Installation complete!"
else
    echo "Unable to install read simulator"
fi