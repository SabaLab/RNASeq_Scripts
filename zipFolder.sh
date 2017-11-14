#!/bin/bash

#############################################################################################
# Spencer Mahaffey
# Oct 2016
# Zip all the files with a given extension in a folder running N in parallel.
# zipFolder.sh /path/to/zip file_extension #_toRunSimultaneously
# 0 - Path to search for files to zip
# 1 - File extension to zip
# 2 - # of files to zip at once
##############################################################################################

PATH1=$1
END=$2
MAX=$3

ZIPLIST=$PATH1"*"$END

COUNT=0

bash_pid=$$

if [ "${#ZIPLIST[@]}" -gt "0" ]
then
        for f in $ZIPLIST
        do
		RUNNINGCOUNT=`ps --no-headers -o pid --ppid=$$ | wc -w`
		RUNNINGCOUNT=$[RUNNINGCOUNT-1]
		echo "$f"
		#echo "running: $RUNNINGCOUNT"
                if [ "$RUNNINGCOUNT" -eq "$MAX" ]
                then
					while [ "$RUNNINGCOUNT" -eq "$MAX" ]
					do
						sleep 60s
						RUNNINGCOUNT=`ps --no-headers -o pid --ppid=$$ | wc -w`
						RUNNINGCOUNT=$[RUNNINGCOUNT-1]
					done
                fi
                #echo "unzipping $f"
                gzip -f6 $f &
        done
	
fi

RUNNINGCOUNT=`ps --no-headers -o pid --ppid=$$ | wc -w`

while [ "$RUNNINGCOUNT" -gt "0" ]
do
        sleep 60s
	RUNNINGCOUNT=`ps --no-headers -o pid --ppid=$$ | wc -w`
	RUNNINGCOUNT=$[RUNNINGCOUNT-1]
done
