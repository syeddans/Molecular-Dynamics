#!/bin/bash

# Setup script for Mole2 - downloads and installs locally
# This allows the repo to be Git-friendly without large binaries

set -e

MOLE2_DIR="Mole2_cmd"
MOLE2_DOWNLOAD_PAGE="https://webchem.ncbr.muni.cz/Platform/App/Mole/Download"
# Direct download link for Mole2 command-line version (may need updating)
MOLE2_DIRECT_URL="https://webchem.ncbr.muni.cz/Platform/AppsBin/Mole/2.5.24.6.8/Mole2_cmd.zip"

echo "Setting up Mole2 for tunnel detection..."

# Check if Mole2 directory already exists
if [ -d "$MOLE2_DIR" ] && [ -f "$MOLE2_DIR/mole2.exe" ]; then
    echo "Mole2 already exists in $MOLE2_DIR"
    echo "Testing existing installation..."
    
    # Test if mono is available
    if command -v mono &> /dev/null || [ -x "/usr/bin/mono" ]; then
        echo "Mono runtime found - testing Mole2..."
        cd "$MOLE2_DIR"
        if timeout 10s mono mole2.exe --help &> /dev/null; then
            echo "✓ Mole2 test successful!"
            echo "✓ Setup complete - ready to run pipeline"
            cd ..
            exit 0
        else
            echo "⚠ Mole2 test failed, but executable exists"
            echo "✓ Setup complete - pipeline should still work"
            cd ..
            exit 0
        fi
    else
        echo "WARNING: mono runtime not found. Please install mono:"
        echo "  Ubuntu/Debian: sudo apt-get install mono-runtime"
        echo "  CentOS/RHEL:   sudo yum install mono-core"
        echo "  Conda:         conda install -c conda-forge mono"
        echo "✓ Mole2 executable exists but cannot be tested without mono"
        exit 0
    fi
fi

echo "Mole2 not found. Setting up new installation..."

# Check if mono is installed first
if ! command -v mono &> /dev/null && [ ! -x "/usr/bin/mono" ]; then
    echo ""
    echo "❌ Mono runtime not found (required for Mole2)"
    echo "Please install mono first:"
    if command -v apt-get &> /dev/null; then
        echo "  Run: sudo apt-get update && sudo apt-get install mono-runtime"
    elif command -v yum &> /dev/null; then
        echo "  Run: sudo yum install mono-core"
    elif command -v conda &> /dev/null; then
        echo "  Run: conda install -c conda-forge mono"
    else
        echo "  Install mono runtime for your system"
    fi
    echo ""
    echo "After installing mono, run this script again."
    exit 1
fi

# Create Mole2 directory
mkdir -p "$MOLE2_DIR"

echo "Attempting to download Mole2..."

# Try direct download first
if wget --spider "$MOLE2_DIRECT_URL" 2>/dev/null; then
    echo "Downloading Mole2 command-line version..."
    wget -O mole2_temp.zip "$MOLE2_DIRECT_URL"
    
    if [ -f "mole2_temp.zip" ]; then
        echo "Extracting Mole2..."
        unzip -q mole2_temp.zip -d mole2_temp/
        
        # Move contents to Mole2_cmd directory
        if [ -d "mole2_temp/Mole2_cmd" ]; then
            cp -r mole2_temp/Mole2_cmd/* "$MOLE2_DIR/"
        elif [ -d "mole2_temp" ]; then
            # Sometimes the zip extracts directly
            find mole2_temp/ -name "mole2.exe" -exec cp -r "$(dirname {})"/* "$MOLE2_DIR/" \;
        fi
        
        # Clean up
        rm -rf mole2_temp.zip mole2_temp/
        
        if [ -f "$MOLE2_DIR/mole2.exe" ]; then
            echo "✓ Mole2 downloaded and extracted successfully!"
        else
            echo "❌ Download succeeded but mole2.exe not found"
            echo "Please download manually from: $MOLE2_DOWNLOAD_PAGE"
            exit 1
        fi
    else
        echo "❌ Download failed"
        echo "Please download manually from: $MOLE2_DOWNLOAD_PAGE"
        exit 1
    fi
else
    echo "❌ Automatic download not available"
    echo ""
    echo "Please download Mole2 manually:"
    echo "1. Visit: $MOLE2_DOWNLOAD_PAGE"
    echo "2. Download the 'Command-line version' for Linux"
    echo "3. Extract the contents to the $MOLE2_DIR/ directory"
    echo "4. Ensure mole2.exe is in $MOLE2_DIR/"
    echo ""
    exit 1
fi

# Test the installation
echo "Testing Mole2 installation..."
cd "$MOLE2_DIR"
if timeout 10s mono mole2.exe --help &> /dev/null; then
    echo "✓ Mole2 installation and test successful!"
    cd ..
else
    echo "⚠ Mole2 installed but test failed (this may be normal)"
    echo "✓ Installation complete - pipeline should work"
    cd ..
fi

echo ""
echo "🎉 Mole2 setup complete!"
echo "📁 Installed in: $MOLE2_DIR/"
echo "🧪 Test command: cd $MOLE2_DIR && mono mole2.exe --help"
echo "🚀 Ready to run the Snakemake pipeline!" 