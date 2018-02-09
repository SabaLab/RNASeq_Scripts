#!/usr/local/bin/python3.5

#############################################################################################
# Spencer Mahaffey
# Jan 2018
# Generic script to run cutadapt on a batch of samples
#
# 1.0.0 - 1/4/18 
# 1.0.1 - 2/9/18 - fixed a bug where unpaired data still required a pair label
#		- fixed a bug parsing sample names when the sample delimiter is not present
##############################################################################################

import sys,subprocess,glob,os

from optparse import OptionParser
from multiprocessing import Process, Queue

def runProcessing(q,cutadaptPath,options,opts):
	while(not q.empty()):
		l=q.get()
		sampleNameWPath=l[0]
		outputFile=l[1]
		sampleName=l[2]
		#compile file list
		if(opts.countRaw):
			runRawCount(sampleNameWPath)
		search=sampleNameWPath+"*"+opts.inSuffix
		if(opts.paired):
			search=sampleNameWPath+"*"+opts.pairPrefix+"1*"+opts.inSuffix
		files=glob.glob(search)
		runCutadapt(cutadaptPath,options,opts,files,sampleName)
		runCountTrimmed(opts,sampleName)

def runRawCount(sampleName):
	print("counting RawReads")
	fileList=glob.glob(sampleName+"*"+opts.inSuffix)
	for f in fileList:
		countFile=f+".count.txt"
		if( not(os.path.isfile(countFile)) or os.stat(countFile).st_size == 0 ):
			countCmd="gzip -dc "+f+" | awk '/@/ {getline; print length($0);getline;getline}' | awk -v sample=\""+f+"\" '{sum+=$1} END {print sample,sum/NR,NR}' > "+countFile
			print("running:"+countCmd)
			completed=subprocess.run(countCmd,shell=True, executable='/bin/bash')
	

def runCutadapt(cutadaptProg,options,opts,files,sampleName):
	print("running cutadapt")
	stdOut=opts.output+"/"+sampleName+".out"
	part1=""
	part2=""
	opart1=" -o "+opts.output+"/"+sampleName+opts.pairPrefix+"1"+".trimmed"+opts.inSuffix
	opart2=" -p "+opts.output+"/"+sampleName+opts.pairPrefix+"2"+".trimmed"+opts.inSuffix
	for f in files:
		if(len(part1)>0):
			part1+=","
			part2+=","
		part1+=f
		if(opts.paired):
			f2=f.replace(opts.pairPrefix+"1",opts.pairPrefix+"2")
			part2+=f2

	fileString=opart1
	if(opts.paired):
		fileString+=opart2
	fileString+=" "+part1+" "
	if(opts.paired):
		fileString=fileString+part2+" "
	cmd=cutAdaptProg+options+fileString+" > "+stdOut+" 2>&1"
	print("running "+cmd)
	completed=subprocess.run(cmd,shell=True, executable='/bin/bash')

def runCountTrimmed(opts,sampleName):
	print("counting trimmedReads")
	files=glob.glob(opts.output+"/"+sampleName+"*.trimmed"+opts.inSuffix)
	for f in files:
		countFile=f+".count.txt"
		countCmd="gzip -dc "+f+" | awk '/@/ {getline; print length($0);getline;getline}' | awk -v sample=\""+f+"\" '{sum+=$1} END {print sample,sum/NR,NR}' > "+countFile
		print("running:"+countCmd)
		completed=subprocess.run(countCmd,shell=True, executable='/bin/bash')

if __name__ == '__main__':
	usage = 'USAGE: %prog inputPath\n'
	parser = OptionParser(usage=usage, version="1.0.1")

	parser.add_option('-c', '--count-raw',  action="store_true", dest='countRaw',default=False,
	                      help='When set script will count/summarize the raw files first')

	parser.add_option('', '--cutadapt-path',  action="store", dest='cutadaptPath', default="",
	                      help='set cutadapt path')
	parser.add_option('', '--cutadapt-q',  action="store", dest='cutadaptQ',type="int", default=20,
	                      help='set cutadapt -q to trim based on quality score')
	parser.add_option('', '--cutadapt-m',  action="store", dest='cutadaptm',type="int", default=20,
	                      help='set cutadapt -m to set the minimum read length.')
	parser.add_option('', '--cutadapt-M',  action="store", dest='cutadaptM',type="int",
	                      help='set cutadapt -M to set the maximum read length.')

	parser.add_option('-a', '--cutadapt-a',  action="store", dest='adapta',type="string", 
	                      help='set cutadapt -a set 3\' adapter for reads or read1 if paired')
	parser.add_option('-A', '--cutadapt-A',  action="store", dest='adaptA',type="string", 
	                      help='set cutadapt -A set 3\' adapter for reads or read2 if paired')
	parser.add_option('-b', '--cutadapt-b',  action="store", dest='adaptb',type="string", 
	                      help='set cutadapt -b set 5\' adapter for reads or read1 if paired')
	parser.add_option('-B', '--cutadapt-B',  action="store", dest='adaptB',type="string", 
	                      help='set cutadapt -B set 5\' adapter for reads or read2 if paired')
	parser.add_option('-g', '--cutadapt-g',  action="store", dest='adaptg',type="string", 
	                      help='set cutadapt -g set 3\' or 5\' adapter for reads or read1 if paired')
	parser.add_option('-G', '--cutadapt-G',  action="store", dest='adaptG',type="string", 
	                      help='set cutadapt -G set 3\' or 5\' adapter for reads or read2 if paired')

	parser.add_option('-P', '--paired',  action="store_true", dest='paired',
	                      help='pass rsem the --paired-end parameter and look for paired end files to pass appropriate paired end files to rsem')
	parser.add_option('-U', '--unpaired',  action="store_false", dest='paired', default=False,
	                      help='pass rsem appropriate unpaired files/parameters')
	parser.add_option('-d', '--delim',  action="store", dest='sampleDelim', default="_L00",
	                      help='A delimiter to detect the end of the sample label, default is _L00 to parse everything before the lane as the sample name.')
	parser.add_option('-p', '',  action="store", dest='maxP',type="int", default=4,
	                      help='set number of processes to run at once.')

	parser.add_option('', '--pair-prefix',  action="store", dest='pairPrefix', default="_R",
	                      help='A prefix before the paired label part of the file name for paired reads, defaults to _R and assumes _R1 - first read and _R2 is second read')
	parser.add_option('-o', '--output',  action="store", dest='output', 
	                      help='The output folder.  Output will go to a folder with the extracted Sample Name in this location.')
	parser.add_option('-i', '--input-suffix',  action="store", dest='inSuffix', default=".fq.gz",
	                      help='The suffix to look for in the inpute files.  ex .fq.gz or .fastq')

	#parser.add_option('', '--combine-across-folders',  action="store_true", dest='combine',default=False,
	#                      help='When set the sample name is parsed based on the -d value and all files across input folders that match are used as input for an RSEM run.')



	(opts, args) = parser.parse_args()
	errors = ''
	if len(args) != 1:
		errors="1 parameter is required in addition to any options provided. Run with -h for help."
	if errors:
		parser.error(errors)
	inPathList=args[0].split(",")
	paired=False
	outPath=""
	options=""
	q = Queue()
	sampleList={}

	cutAdaptProg="cutadapt "

	if(opts.cutadaptPath):
		cutAdaptProg=opts.cutadaptPath
		if(not cutAdaptProg[-1]=="/"):
			cutAdaptProg+="/"
		cutAdaptProg+="cutadapt "

	if(opts.cutadaptQ):
		options+=" -q "+str(opts.cutadaptQ)
	if(opts.cutadaptm):
		options+=" -m "+str(opts.cutadaptm)
	if(opts.cutadaptM):
		options+=" -M "+str(opts.cutadaptM)
	if(opts.adapta):
		options+=" -a "+opts.adapta
	if(opts.adaptA):
		options+=" -A "+opts.adaptA
	if(opts.adaptb):
		options+=" -b "+opts.adaptb
	if(opts.adaptB):
		options+=" -B "+opts.adaptB
	if(opts.adaptg):
		options+=" -g "+opts.adaptg
	if(opts.adaptG):
		options+=" -G "+opts.adaptG

	if(opts.paired):
		paired=True
	if(opts.output):
		outPath=opts.output+"/"

	for inPath in inPathList:
		print("processing:"+inPath)
		search="*"+opts.inSuffix
		if(opts.paired):
			search="*"+opts.pairPrefix+"1*"+opts.inSuffix
		print("Search:"+inPath+"/"+search)
		fileList=glob.glob(inPath+"/"+search)
		for f in fileList:
			print("file"+f)
			start=f.rfind("/")+1
			end=f.rfind(opts.sampleDelim)
			#send=f.find("_",start)
			if end<0:
				end=f.rfind(opts.inSuffix)
				if(end<0):
					sampleName=f[start:]
				else:
					sampleName=f[start:end]
			else:
				sampleName=f[start:end]
			outFile=outPath+sampleName+".trimmed.fq.gz"
			print("adding:"+inPath+"/"+sampleName)
			q.put([inPath+"/"+sampleName,outFile,sampleName])
	#run rawCounting/cutadapt in parallel. for opts.p processes at a time.
	pList=[]
	for i in range(0,(opts.maxP)):
		pList.append(Process(target=runProcessing, args=(q,cutAdaptProg,options,opts)))
		pList[i].start()

	for i in range(0,(opts.maxP)):
		pList[i].join()

	sumCmd="cat "+inPath+"/*.count.txt >> "+inPath+"/countSummary.txt"
	print("sumarize raw counts")
	print("running:"+sumCmd)
	completed=subprocess.run(sumCmd,shell=True, executable='/bin/bash')

	#run summary of trimmedReads
	sumCmd="cat "+outPath+"/*.count.txt >> "+outPath+"/trimmedSummary.txt"
	print("sumarize trimmed counts")
	print("running:"+sumCmd)
	completed=subprocess.run(sumCmd,shell=True, executable='/bin/bash')
