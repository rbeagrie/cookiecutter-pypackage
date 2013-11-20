#!/bin/bash

for f in "$@"
do 
    READS=$READS$f
    READS=$READS','
done

echo $READS | sed s/,$//
