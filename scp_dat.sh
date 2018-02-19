#!/bin/bash
if [ $# -lt 1 ]; then
    echo "Usage : scp_dat.sh dat.id"
    exit;
fi

dat_file=$1

scp felisma@137.122.235.201:/oisb/data/*/$dat_file ./


