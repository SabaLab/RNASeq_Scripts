#!/usr/local/bin/python3.5

#############################################################################################
# Spencer Mahaffey
# Oct 2017
# Generic script to run RSEM on a batch of samples with one index or strain specific indices.
#
# 1.1 - 12/22/17 - Add support to combine samples across batches
##############################################################################################

import sys,subprocess,glob,os

from optparse import OptionParser

usage = 'USAGE: %prog inputPath indexPath indexFile numProcessors\n'
parser = OptionParser(usage=usage, version="1.1.0")

parser.add_option('-a', '--input-sample-dir',  action="store_true", dest='inputSampleDir',default=False,
                      help='When set input is assumed to be in /path/sample/unmapped.end1.fq etc.  default assumes /path/sample.fq')
parser.add_option('', '--rsem-time',  action="store_true", dest='rsemTime', default=False,
                      help='pass rsem the --time parameter')
parser.add_option('', '--rsem-seedLen', action="store", dest='rsemSeedLen',type='int', default=0,
                      help='pass rsem the --seed-length parameter and value')
parser.add_option('', '--rsem-seed', action="store", dest='rsemSeed',type='int', default=0,
                      help='pass rsem the --seed parameter and value')
parser.add_option('', '--rsem-bowtie2',  action="store_true", dest='rsemBowtie2', default=False,
                      help='pass rsem the --bowtie2 parameters')
parser.add_option('', '--rsem-noBam', action="store_true", dest='rsemNoBam', default=False,
                      help='pass rsem the --no-bam-ouput parameter')
parser.add_option('', '--rsem-fwProb', action="store", dest='rsemFwProb',
                      help='pass rsem the --forward-prob parameter and value')
parser.add_option('-P', '--paired',  action="store_true", dest='paired',
                      help='pass rsem the --paired-end parameter and look for paired end files to pass appropriate paired end files to rsem')
parser.add_option('-U', '--unpaired',  action="store_false", dest='paired', default=False,
                      help='pass rsem appropriate unpaired files/parameters')
parser.add_option('-d', '--delim',  action="store", dest='sampleDelim', default="_L00",
                      help='A delimiter to detect the end of the sample label, default is _L00 to parse everything before the lane as the sample name.')
parser.add_option('', '--pair-prefix',  action="store", dest='pairPrefix', default="_R",
                      help='A prefix before the paired label part of the file name for paired reads, defaults to _R and assumes _R1 - first read and _R2 is second read')
parser.add_option('-o', '--output',  action="store", dest='output', 
                      help='The output folder.  Output will go to a folder with the extracted Sample Name in this location.')
parser.add_option('-i', '--input-suffix',  action="store", dest='inSuffix', 
                      help='The suffix to look for in the inpute files.  ex .fq.gz or .fastq')
parser.add_option('-s', '--index-ssg',  action="store_true", dest='ssgIndex',default=False,
                      help='When set the begining of the filename is expected to denote strain.  The strain will be parsed and used to align to a strain specific genome.')
parser.add_option('', '--combine-across-folders',  action="store_true", dest='combine',default=False,
                      help='When set the sample name is parsed based on the -d value and all files across input folders that match are used as input for an RSEM run.')
parser.add_option('', '--output-pbs',  action="store_true", dest='pbs',default=False,
                      help='When set instead of calling RSEM directly output goes to a pbs file')
parser.add_option('', '--pbs',  action="store", dest='pbsFile',default="",
                      help='PBS output file name.')


(opts, args) = parser.parse_args()
errors = ''
if len(args) != 4:
	errors="4 parameters are required in addition to any options provided. Run with -h for help."
if errors:
	parser.error(errors)
inPathList=args[0].split(",")
indexPath=args[1]
indexSuffix=args[2]
maxP=args[3]
paired=False
outPath=""

options=""

sampleList={}
PBS=""

if(opts.rsemTime):
	options+=" --time"
if(opts.rsemSeed>0):
	options+=" --seed "+str(opts.rsemSeed)
if(opts.rsemSeedLen>0):
	options+=" --seed-length "+str(opts.rsemSeedLen)
if(opts.rsemBowtie2):
	options+=" --bowtie2"
if(opts.rsemNoBam):
	options+=" --no-bam-output"
if(opts.rsemFwProb):
	options+=" --forward-prob "+opts.rsemFwProb
if(opts.paired):
	options+=" --paired-end"
	paired=True
if(opts.output):
	outPath=opts.output+"/"

if(opts.pbs):
	PBS=open(opts.pbsFile,"w")
	PBS.write("#!/bin/sh\n")
	PBS.write("#PBS -l select=1:ncpus=16:nmics=1:mem=112gb\n")
	PBS.write("#PBS -l walltime=700:00:00\n")
	PBS.write("#PBS -o rsemb\n")
	PBS.write("#PBS -N rsemb\n\n")

	PBS.write("mkdir /state/partition1/mahaffey")
	options+=" --temporary-folder /state/partition1/mahaffey/tmp"

for inPath in inPathList:
	search=inPath+"/*"+opts.pairPrefix+"1*"+opts.inSuffix
	if(not paired):
		search=inPath+"/*"+opts.inSuffix
	if(opts.inputSampleDir):
		search=inPath+"/*/"
	print(search)
	fileList=glob.glob(search)
	for f in fileList:
		start=f.rfind("/")+1
		end=f.rfind(opts.sampleDelim)
		if(opts.inputSampleDir):
			start=f.rfind("/",0,f.rfind("/")-2)+1
		if end==-1:
			sampleName=f[start:]
		else:
			sampleName=f[start:end]
		outSample=outPath+sampleName
		if(not sampleName in sampleList):
			sampleList[sampleName]=[]
			index=indexPath+"/"+indexSuffix
			if(opts.ssgIndex):
				firstDelim="_"
				underscorePos=f.find("_",start)
				hyphenPos=f.find("-",start)
				if(hyphenPos>-1 and underscorePos>-1 and hyphenPos<underscorePos):
					firstDelim="-"
				elif(hyphenPos>-1 and underscorePos==-1):
					firstDelim="-"
				send=f.find(firstDelim,start)
				strain=f[start:send]
				if strain=='Dark':
					strain="DA"
				elif strain=='ACI':
					strain="ACI_EurMcwi"
				elif strain=='Cop':
					strain="BN"
				index=indexPath+"/"+strain+indexSuffix
		    
			if(not os.path.exists(outSample+".stat")):
				unmappedR1=""
				unmappedR2=""
				ucount=0
				for inPath2 in inPathList:
					if(opts.combine or (not opts.combine and inPath2==inPath)):
						found=False
						search2=inPath2+"/"+sampleName+"*"+opts.pairPrefix+"1*"+opts.inSuffix
						if(not paired):
							search2=inPath2+"/"+sampleName+"*"+opts.inSuffix
						if(opts.inputSampleDir):
							search2=inPath2+"/"+sampleName+"*/"
						flist=glob.glob(search2)
						for f2 in flist:
							found=True
							if(ucount>0):
								unmappedR1=unmappedR1+","
								if(paired):
									unmappedR2=unmappedR2+","
							if(opts.inputSampleDir):
								unmappedR1=unmappedR1+f2+"/unmapped.end1"+opts.inSuffix
								if(paired):
									unmappedR2=unmappedR2+f2+"/unmapped.end2"+opts.inSuffix
							else:
								unmappedR1=unmappedR1+f2
								if(paired):
									unmappedR2=unmappedR2+f2.replace(opts.pairPrefix+"1",opts.pairPrefix+"2")
							ucount+=1
						if(not found):
							print("WARNING: "+sampleName+" not found in "+inPath2)
				arg=str("rsem-calculate-expression -p "+maxP+options+" "+unmappedR1)
				if(paired):
					arg=arg+str(" "+unmappedR2)
				arg=arg+str(" "+index+" "+outSample+" > "+sampleName+".rsem.out 2>&1")
				if(opts.pbs):
					PBS.write("echo 'running: "+sampleName+"'\n")
					PBS.write(arg+"\n")
					PBS.write("rm -r /state/partition1/mahaffey/tmp\n")
				else:
					print("running: "+sampleName)
					print(arg)
					print("\n")
					completed=subprocess.run(arg,shell=True)
if(opts.pbs):
	PBS.write("rm -r /state/partition1/mahaffey\n")
	PBS.close()
