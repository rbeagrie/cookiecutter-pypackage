#!/bin/bash

db=`ls *.chrom.sizes | sed 's/.chrom.sizes//'`
echo $db

wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes
chmod +x fetchChromSizes
./fetchChromSizes $db > $db.chrom.sizes

