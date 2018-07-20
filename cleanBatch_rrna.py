#!/usr/bin/env python

#############################################################################################
# Spencer Mahaffey
# April 2018
# Generic script to run bowtie2 alignment to rRNA on a batch of samples.
#
# 1.0 - 4/3/18 - Initial version
# 1.0.1 - 7/20/18 - Fixed Sytax error, add to help
##############################################################################################

import sys,subprocess,glob,os

from optparse import OptionParser

usage = 'USAGE: %prog inputPath indexPrefix procNum\n'
parser = OptionParser(usage=usage, version="1.0.1")
parser.add_option('-p', '',  action="store", dest='maxProc',type='int',default=1,
                      help='set the number of processors to run aligner with')

parser.add_option('-P', '--paired',  action="store_true", dest='paired',
                      help='pass rsem the --paired-end parameter and look for paired end files to pass appropriate paired end files to rsem')

parser.add_option('-U', '--unpaired',  action="store_false", dest='paired', default=False,
                      help='pass rsem appropriate unpaired files/parameters')

parser.add_option('-d', '--delim',  action="store", dest='sampleDelim', default="_L00",
                      help='A delimiter to detect the end of the sample label, default is _L00 to parse everything before the lane as the sample name.')
parser.add_option('', '--pair-prefix',  action="store", dest='pairPrefix', default="_R",
                      help='A prefix before the paired label part of the file name for paired reads, defaults to _R and assumes _R1 - first read and _R2 is second read')

parser.add_option('-i', '--input-suffix',  action="store", dest='inSuffix', default=".fq.gz",
                      help='The suffix to look for in the inpute files.  ex .fq.gz or .fastq')

parser.add_option('-o', '--output',  action="store", dest='output', 
                      help='The output folder.  Output will go to a folder with the extracted Sample Name in this location.')

parser.add_option('-x', '--index',  action="store",type="string", dest='index',help='index directory')
parser.add_option('-F', '--unmapped-fastq',  action="store_false", dest='toFastQ', default=False,
                      help='convert unmapped.bam to fastq files')
parser.add_option('', '--bt2-input',  action="store", dest='bt2Input', default="q",
                      help='Bowtie2 Input file format (default -q)')
parser.add_option('', '--bt2-k',  action="store", dest='bt2K', default=2,
                      help='Bowtie2 option -k ## (default 2)')
parser.add_option('', '--bt2-sensitivity',  action="store", dest='bt2Sensitivity', default="very-sensitive",
                      help='Bowtie2 option (default --very-sensitive)')
parser.add_option('', '--bt2-gbar',  action="store", dest='bt2Gbar', default=4,
                      help='Bowtie2 option --gbar(#) default(4)')
parser.add_option('', '--bt2-mp',  action="store", dest='bt2Mp', default="10,4",
                      help='Bowtie2 --mp #,# default(10,4)')
parser.add_option('', '--bt2-np',  action="store", dest='bt2Np', default=1,
                      help='Bowtie2 --np # default(1)')
parser.add_option('', '--bt2-stranded',  action="store", dest='bt2Stranded', default="fr",
                      help='Bowtie2 stranded option default(--fr)')


(opts, args) = parser.parse_args()
errors = ''
if len(args) != 1:
    errors="1 argument is required in addition to any options provided. Run with -h for help."
if errors:
    parser.error(errors)

inPath=args[0]
options=""

if(opts.bt2Input):
    options+="-"+opts.bt2Input+" "
if(opts.bt2K):
    options+="-k "+str(opts.bt2K)+" "
if(opts.bt2Sensitivity):
    options+="--"+opts.bt2Sensitivity+" "
if(opts.bt2Gbar):
    options+="--gbar "+str(opts.bt2Gbar)+" "
if(opts.bt2Mp):
    options+="--mp "+opts.bt2Mp+" "
if(opts.bt2Np):
    options+="--np "+str(opts.bt2Np)+" "
if(opts.bt2Stranded):
    options+="--"+opts.bt2Stranded+" "
if(opts.maxProc>1):
    options+="-p "+str(opts.maxProc)+" "

search="*"
if(opts.paired):
    search="*"+opts.pairPrefix+"1*"

fileList=glob.glob(inPath+"/"+search+opts.inSuffix)

for f in fileList:
    sampleName=""
    start=f.rfind("/")+1
    end=f.rfind(opts.sampleDelim)-1
    if(inPath):
        start=f.rfind("/",0,f.rfind("/")-1)  
    if end==-1:
        sampleName=f[start:]
    else:
        sampleName=f[start:end]
    unmappedBAM=opts.output+"/"+sampleName+"/unmapped.bam"
    alignedBAM=opts.output+"/"+sampleName+"/aligned.bam"
    if(not os.path.exists(opts.output+"/"+sampleName+"/aligned.bam")):
        if(not os.path.exists(opts.output+"/"+sampleName)):
            os.makedirs(opts.output+"/"+sampleName)
        f2=f.replace(opts.pairPrefix+"1",opts.pairPrefix+"2")
        arg=str("bowtie2 "+options+"-x "+opts.index+" -1 "+f)
        if(opts.paired):
            arg+=" -2 "+f2
        arg+=" 2> "+opts.output+"/"+sampleName+"/alignSummary.txt | tee >(samtools view -f 0x4 -bS1 - > "+unmappedBAM+" 2> /dev/null) >(samtools view -F 0x4 -bS1 - > "+alignedBAM+" 2> /dev/null) | grep errors "
        print("running: "+sampleName)
        print("w command: "+arg)
        completed=subprocess.run(arg,shell=True, executable='/bin/bash')

        if(opts.toFastQ):
            convertBam="samtools fastq "+unmappedBAM+" | gzip > "+opts.output+"/"+sampleName+"/"+unmapped.end1.fq.gz
            if(opts.paired):
                convertBam="/usr/local/scripts/bamToFastQ_py.paired.sh"+opts.output+"/"+sampleName+" "+unmappedBAM
            completed=subprocess.run(convertBam,shell=True, executable='/bin/bash')

