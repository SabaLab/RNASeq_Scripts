#!/bin/bash

#############################################################################################
# Spencer Mahaffey
# Apr 2016
# Automate converting bam files to stranded bigWig files for UCSC genome browser.
##############################################################################################

if [ "$#" -lt 4 ]; then
        echo "Usage: bamToBigWig.Stranded.sh /path/to/DirWBam bamFileName /path/to/BigWigOutput chromSizesFile"
        echo ""
	echo "Convert stranded bam alignments file into 2 bigWig files. One for the plus strand"
        echo "one minus strand. Creates additional temporary files in the working directory with the"
        echo "bam file. Note: If your chrom.sizes start with chr and intermediate bedGraph files"
        echo "don't you need use bamToBigWig.Stranded.wChr.sh as this will insert 'chr' on each line."
	echo ""
        echo "1- Working Directory"
        echo "2- Bam file"
        echo "3- BigWig Output Path"
        echo "4- Chrom.sizes file"
        echo ""
        echo ""
        exit 0
fi


INPUTDIR=$1
INPUTFILE=$2
FINALOUT=$3
CHROMSIZES=$4


len=${#INPUTFILE}-4
INPUTFILEPART=${INPUTFILE:0:$len}


samtools view -b -f 0X40 $INPUTDIR$INPUTFILE | tee >(bedtools genomecov -bg -split -strand + -ibam - > $INPUTFILEPART.first.plus.bg) >(bedtools genomecov -bg -split -strand - -ibam - > $INPUTFILEPART.first.minus.bg) | grep errors &
S1=$!
samtools view -b -f 0X80 $INPUTDIR$INPUTFILE | tee >(bedtools genomecov -bg -split -strand + -ibam - > $INPUTFILEPART.second.plus.bg) >(bedtools genomecov -bg -split -strand - -ibam - > $INPUTFILEPART.second.minus.bg) | grep errors &
S2=$!

wait $S1
wait $S2

bedtools unionbedg -i $INPUTFILEPART.first.minus.bg $INPUTFILEPART.second.plus.bg | awk '{ print $1 "\t" $2 "\t" $3 "\t" ($4+$5) }' - > $INPUTFILEPART.plus.combined.bg &
BG1=$!
bedtools unionbedg -i $INPUTFILEPART.first.plus.bg $INPUTFILEPART.second.minus.bg | awk '{ print $1 "\t" $2 "\t" $3 "\t" ($4+$5) }' - > $INPUTFILEPART.minus.combined.bg &
BG2=$!
wait $BG1
wait $BG2

awk '{if($1~"ERCC") print $0}' $INPUTFILEPART.plus.combined.bg > $INPUTFILEPART.plus.combined.ercc.bg &
awk '{if($1~"ERCC") print $0}' $INPUTFILEPART.minus.combined.bg > $INPUTFILEPART.minus.combined.ercc.bg &
awk '{if($1!~"ERCC") print $0}' $INPUTFILEPART.plus.combined.bg > $INPUTFILEPART.plus.combined.no_ercc.bg &
A3=$!
awk '{if($1!~"ERCC") print $0}' $INPUTFILEPART.minus.combined.bg > $INPUTFILEPART.minus.combined.no_ercc.bg &
A4=$!

wait $A3
wait $A4

#sed -i -e 's/^/chr/g' $INPUTFILEPART.plus.combined.no_ercc.bg &
#SEDPlus=$!
#sed -i -e 's/^/chr/g' $INPUTFILEPART.minus.combined.no_ercc.bg &
#SEDMinus=$!

wait $SEDPlus
wait $SEDMinus

bedGraphToBigWig $INPUTFILEPART.plus.combined.no_ercc.bg  $CHROMSIZES  $FINALOUT$INPUTFILEPART.plus.bw &
BGW1=$!
bedGraphToBigWig $INPUTFILEPART.minus.combined.no_ercc.bg  $CHROMSIZES  $FINALOUT$INPUTFILEPART.minus.bw &
BGW2=$!

wait $BGW1
wait $BGW2
