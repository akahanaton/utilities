#!/bin/bash
if [ $# -lt 1 ]; then
    echo "Usage : xtandem_parser.sh xml_file"
    exit;
fi

xml_file=$1

R --slave --no-restore --no-save --no-readline $xml_file < ~/myProgram/xtandem_parser.R
