#!/bin/bash 
COUNTER=0
while [ $COUNTER -le 10 ]; do
    echo The counter is $COUNTER
    let COUNTER=COUNTER+1 
done
