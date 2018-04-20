#!/usr/bin/env python

#############################################################################################
# Spencer Mahaffey
# Mar 2018
# Generic script to run stringtie on a batch of samples.
#
# 1.0 - 3/29/18 - Initial version
##############################################################################################

import sys,subprocess,glob,os

from optparse import OptionParser

usage = 'USAGE: %prog options inputPath\n'
parser = OptionParser(usage=usage, version="1.0.0")

parser.add_option('-a', '--input-sample-dir',  action="store_true", dest='inputSampleDir',default=False,
                      help='When set input is assumed to be in /path/sample/unmapped.end1.fq etc.  default assumes /path/sample.fq')

parser.add_option('-p', '',  action="store", dest='maxProc',type='int',default=1,
                      help='set the number of processors to run stringtie with')

parser.add_option('-o', '--output',  action="store", dest='output', 
                      help='The output folder.  Output will go to a folder with the extracted Sample Name in this location.')

parser.add_option('', '--st-v',  action="store_true", dest='stV', default=False,
                      help='pass stringtie the -v parameter')

parser.add_option('-g', '--st-G', action="store", dest='gtf', default="",
                      help='pass stringtie the -G parameter for guided reconstruction with GTF')

parser.add_option('-l', '--st-l', action="store", dest='label', default="",
                      help='pass stringtie the -l parameter and value - label for output')

parser.add_option('-s', '',  action="store_true", dest='sort', default=False,
                      help='sort bam files')

parser.add_option('-m', '', action="store_true", dest='merge', default=False,
                      help='merge bam files')

parser.add_option('', '--s-p',  action="store", dest='sortProc',type='int',default=1,
                      help='set the number of processors to run sort')

parser.add_option('', '--s-M',  action="store", dest='sortMem',default="1G",
                      help='set the memory per thread to run sort with')

#parser.add_option('', '--parallel',  action="store_true", dest='maxProc',default=false,
#                      help='set the number of processors to run stringtie with')

parser.add_option('', '--output-pbs',  action="store_true", dest='pbs',default=False,
                      help='When set instead of calling Stringtie directly output goes to a pbs file')

parser.add_option('', '--pbs',  action="store", dest='pbsFile',default="",
                      help='PBS output file name.')

parser.add_option('-d', '', action="store", dest='delim', default="_L",
                      help='sampleName delimiter to parse sampleNames to match across folders for merging.')

(opts, args) = parser.parse_args()
errors = ''
if len(args) != 1:
  errors="1 parameter is required in addition to any options provided. Run with -h for help."
if errors:
  parser.error(errors)

inPathList=args[0].split(",")
outPath=""
options=""
sampleList={}
PBS=""

if(opts.gtf):
  options+="-G "+opts.gtf+" "
if(opts.stV):
  options+="-v "


if(opts.pbs):
  PBS=open(opts.pbsFile,"w")
  PBS.write("#!/bin/sh\n")
  PBS.write("#PBS -l select=1:ncpus=16:nmics=1:mem=112gb\n")
  PBS.write("#PBS -l walltime=700:00:00\n")
  PBS.write("#PBS -o strtie\n")
  PBS.write("#PBS -N strtie\n\n")

#identify samples and build a list of files for each sample
for inPath in inPathList:
  search=inPath+"/*.bam"
  if(opts.inputSampleDir):
    search=inPath+"/*/aligned.bam"
  print("search:"+search)
  fileList=glob.glob(search)
  for f in fileList:
    start=f.rfind("/")+1
    end=-4 #remove .bam
    if(opts.inputSampleDir):
      start=f.rfind("/",0,f.rfind("/")-2)+1
      end=f.rfind("/",start)+1
    if end==-1:
      sampleName=f[start:]
    else:
      sampleName=f[start:end]
    if(sampleName.find(opts.delim)):
      sampleName=sampleName[0:sampleName.find(opts.delim)]
    if(not sampleName in sampleList):
      sampleList[sampleName]=[]
      sampleList[sampleName].append(sampleName)
    sampleList[sampleName].append(f)


for sampleFileList in sampleList:
  #find occurances of sample
  toMerge=""
  inputFile=""
  sampleName=""
  sOptions=options
  for sampleFile in sampleList[sampleFileList]:
    if(sampleName==""):
      sampleName=sampleFile
    else:
      #sort
      if(opts.sort):
        sopt=""
        if(opts.sortProc>1):
          sopt="-@ "+str(opts.sortProc)+" "
        if(opts.sortMem):
          sopt+="-m "+opts.sortMem+" "
        sortbam=sampleFile.replace(".bam",".sorted.bam")
        arg="samtools sort "+sopt+" -o "+sortbam+" "+sampleFile
        if(opts.pbs):
          PBS.write("echo 'sorting: "+sampleFile+"'\n")
          PBS.write(arg+"\n")
        else:
          print("sorting: "+sampleFile)
          print(arg)
          print("\n")
          completed=subprocess.run(arg,shell=True)
        if(toMerge!=""):
          toMerge+=" "
          inputFile+=","
        toMerge+=sortbam
        inputFile+=sortbam
      else:
        if(toMerge!=""):
          toMerge+=" "
          inputFile+=","
        toMerge+=sampleFile
        inputFile+=sampleFile
  #merge sampleFiles
  if(opts.merge):
    inputFile=opts.output+"/"+sampleName+".merged.bam"
    arg="samtools merge "+opts.output+"/"+sampleName+".merged.bam "+toMerge
    if(opts.pbs):
      PBS.write("echo 'merging: "+sample+"'\n")
      PBS.write(arg+"\n")
    else:
      print("merging: "+sample)
      print(arg)
      print("\n")
      completed=subprocess.run(arg,shell=True)

  #stringtie
  fileList=inputFile.split(",")
  for f in fileList:
    name=sampleName
    if(len(fileList)>1):
      if(opts.inputSampleDir):
        start=f.rfind("/",0,f.rfind("/"))+1
        end=f.rfind("/")
        name=f[start:end]
      else:
        start=f.rfind("/")+1
        end=f.rfind(".bam")
        end2=f.rfind(".sorted.bam")
        low=end
        if(end2>-1):
          low=end2
        name=f[start:low]
    if(opts.output):
      sOptions+= "-o "+opts.output+"/"+name+".gtf "
    if(opts.label):
      sOptions+="-l "+name+opts.label+" "
    arg="stringtie "+sOptions+" "+f+" > "+opts.output+"/"+name+".out 2>&1"

    if(opts.pbs):
      PBS.write("echo 'running: "+sampleName+"'\n")
      PBS.write(arg+"\n")
    else:
      print("running: "+sampleName)
      print(arg)
      print("\n")
      completed=subprocess.run(arg,shell=True)
