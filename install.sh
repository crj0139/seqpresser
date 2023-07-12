#!/bin/bash

# Clone the repository
git clone https://github.com/username/repository.git
cd repository

# Create conda environment from environment.yml
conda env create -f environment.yml

# Activate the newly created environment
conda activate seqpresser

# Set up any dependencies
# Example: Install Python dependencies using pip
pip install -r requirements.txt

# Perform additional setup steps
# Example: Clone seqrequester repository and build
git clone https://github.com/marbl/seqrequester.git
cd seqrequester/src
make

# Return to the repository directory
cd ../../

# Perform any remaining setup steps
# Example: Build and install the project
python setup.py build
python setup.py install

# Optional: Clean up
cd ..
rm -rf repository

echo "Installation complete!"
