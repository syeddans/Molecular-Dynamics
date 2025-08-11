#!/usr/bin/env python3
"""
Fix Mole2 XML by properly quoting attribute values.
"""

import sys
import re

def fix_xml_attributes(xml_content):
    """Fix XML attributes by properly quoting them."""
    # Fix SequenceNumber attributes
    xml_content = re.sub(r'SequenceNumber=(\d+)', r'SequenceNumber="\1"', xml_content)
    # Fix Chain attributes
    xml_content = re.sub(r'Chain=([A-Z])', r'Chain="\1"', xml_content)
    return xml_content

def main():
    if len(sys.argv) != 3:
        print("Usage: python fix_mole2_xml.py <input_xml> <output_xml>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    with open(input_file, 'r') as f:
        content = f.read()
    
    fixed_content = fix_xml_attributes(content)
    
    with open(output_file, 'w') as f:
        f.write(fixed_content)
    
    print(f"Fixed XML saved to {output_file}")

if __name__ == "__main__":
    main() 