#!/bin/env bash
set –euo pipefail
if [ $# -lt 1 ]; then
    echo "Usage : ass_result.sh txt_file sep_tag,default: \t"
    exit;
fi

txt_file=$1
if [ ! -z "$2" ]; then
    sep_tag=$2
else
    sep_tag="\t"
fi

head -2 $txt_file | awk -F "$sep_tag" '{for(i=1;i<=NF;i++) {print i,$i}}'

