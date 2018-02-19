#!/bin/sh
if [ $# -lt 1 ]; then
    echo ""
    echo "\tUsage : $0 infile outfile"
    echo "\tcheck coding type in vi: :set fileencoding"
    echo ""
    exit;
fi

iconv -t utf-8 -f cp936 -c $1 > tmp.file
mv tmp.file $2

