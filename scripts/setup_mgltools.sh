#!/bin/bash

# Setup script for MGLTools - downloads and installs locally
# This allows the repo to be Git-friendly without large binaries

set -e

MGLTOOLS_DIR="mgltools_x86_64Linux2_1.5.6"
MGLTOOLS_URL="https://ccsb.scripps.edu/mgltools/download/495/"

echo "Setting up MGLTools for AutoDock..."

# Download MGLTools if not already present
if [ ! -f "mgltools_x86_64Linux2_1.5.6.tar.gz" ]; then
    echo "Downloading MGLTools..."
    wget --max-redirect=2 --trust-server-names -O mgltools_x86_64Linux2_1.5.6.tar.gz "$MGLTOOLS_URL"
fi

echo "Extracting MGLTools..."
tar -xzf mgltools_x86_64Linux2_1.5.6.tar.gz


# Run the install script that comes with MGLTools
echo "Running MGLTools install.sh..."
pwd
cd "$MGLTOOLS_DIR"
pwd
if [ -f "install.sh" ]; then
    bash install.sh
else
    echo "ERROR: install.sh not found in MGLTools directory"
    exit 1
fi
cd ..

# Check if installation was successful
if [ -f "$MGLTOOLS_DIR/bin/pythonsh" ]; then
    echo "MGLTools installation complete!"
    echo "Python 2.5 interpreter available at: $MGLTOOLS_DIR/bin/pythonsh"
    
    # Test the installation
    echo "Testing MGLTools installation..."
    if "$MGLTOOLS_DIR/bin/pythonsh" -c "print 'MGLTools Python 2.5 working'"; then
        echo "MGLTools test successful!"
    else
        echo "WARNING: MGLTools test failed, but binaries exist"
    fi
else
    echo "ERROR: MGLTools installation failed - pythonsh not found"
    exit 1
fi

# Clean up download
rm -f mgltools_x86_64Linux2_1.5.6.tar.gz

echo "Setup complete. You can now run the Snakemake pipeline." 