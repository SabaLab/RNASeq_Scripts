#!/bin/bash

#############################################################################################
# Spencer Mahaffey
# July 2016
# Line count gzipped files in a folder with a specific file suffix and leaving the zipped files
# in place.
##############################################################################################

CPATH=$1
SUFFIX=$2
LIST=$CPATH"*"$SUFFIX

for f in $LIST
do
        echo "running $f"
	echo $f >> $CPATH"count.txt"
        gunzip -c $f | wc -l >> $CPATH"count.txt"
done

