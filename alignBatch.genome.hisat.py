#!/usr/local/bin/python3.5

#############################################################################################
# Spencer Mahaffey
# Oct 2017
# Generic script to align batches of RNASeq to generic or strain specific genomes with hisat2
##############################################################################################

import sys,subprocess,glob,os

from optparse import OptionParser

usage = 'USAGE: %prog [options] inputDir\n'
parser = OptionParser(usage=usage, version="1.0.2")

parser.add_option('-a', '--input-sample-dir',  action="store_true", dest='inputSampleDir',default=False,
                      help='When set input is assumed to be in /path/sample/unmapped.end1.fq etc.  default assumes /path/sample.fq')

parser.add_option('-o', '--output',  action="store",type="string", dest='outDir',
        help='output directory')

parser.add_option('-x', '--index-dir',  action="store",type="string", dest='indexDir',
        help='index directory')

parser.add_option('', '--index-prefix',  action="store",type="string", dest='indexPref', default="",
        help='index files prefix')

parser.add_option('', '--input-suffix',  action="store",type="string", dest='inSuffix',default=".fq.gz",
        help='input suffix ex default .fq.gz')

parser.add_option('-s', '--strain-specific',  action="store_true", dest='strainSpecific',default=False,
                      help='Enables extraction of the strain from the beginging of the file name and appending to the index-suffix to generate the index file for alignment')

parser.add_option('', '--index-suffix',  action="store",type="string", dest='indexSuffix',
                      help='strain specific index files suffix')

#parser.add_option('','--groupSamples', dest='group', help='Allows grouping files from the same sample in the alignment -- NOT IMPLEMENTED YET')

parser.add_option('-P', '--paired',  action="store_true", dest='paired',
                      help='pass rsem the --paired-end parameter and look for paired end files to pass appropriate paired end files to rsem')

parser.add_option('-U', '--unpaired',  action="store_false", dest='paired', default=False,
                      help='pass rsem appropriate unpaired files/parameters')

parser.add_option('-p', '',  action="store", dest='proc', type='int',
                      help='number of processors to use with hisat2')

parser.add_option('', '--hisat2-k',action="store",type="int", dest='hisatK',default=20,
                      help='hasat2 -k option - number of alignments before stopping')

parser.add_option('', '--hisat2-strandness',action="store",type="string", dest='hisatStrandness',
                      help='hasat2 -rna-strandness option - read strandness')

parser.add_option('', '--hisat2-reorder',action="store_true", dest='hisatReorder',default=False,
                      help='hasat2 -reorder option - reorder output reads')

parser.add_option('', '--hisat2-summary',action="store_true", dest='hisatSummary',default=True,
                      help='hasat2 --new-summary --summary-file options - output summary file')

parser.add_option('', '--hisat2-splice',action="store_true", dest='hisatSplice',default=True,
                      help='hasat2 --novel-splicesite-outfile option - output splice junction file')

parser.add_option('', '--hisat2-input',action="store",type="string", dest='hisatIn',default="q",
                      help='hasat2 input file type options - passes letter as the option to specify file type default is q for fastq.')

parser.add_option('-d', '--delim',  action="store", dest='sampleDelim', default="_L00",
                      help='A delimiter to detect the end of the sample label, default is _L00 to parse everything before the lane as the sample name.')

parser.add_option('', '--pair-prefix',  action="store", dest='pairPrefix', default="_R",
                      help='A prefix before the paired label part of the file name for paired reads, defaults to _R and assumes _R1 - first read and _R2 is second read')




(opts, args) = parser.parse_args()
errors = ''
if len(args) != 1:
        errors="1 parameter is required in addition to any options provided. Run with -h for help."
if errors:
        parser.error(errors)

inPath=args[0]
outPath=inPath
if(opts.outDir):
        outPath=opts.outDir
indexPath=inPath
if(opts.indexDir):
        indexPath=opts.indexDir
maxP=4
if(opts.proc):
        maxP=opts.proc
search="*"+opts.inSuffix
if(opts.paired):
        search="*"+opts.pairPrefix+"1*"+opts.inSuffix
if(opts.inputSampleDir):
        search="*/"
fileList=glob.glob(inPath+"/"+search)
for f in fileList:
        start=f.rfind("/")+1
        end=f.rfind(opts.sampleDelim)-1
        #send=f.find("_",start)
        if(opts.inputSampleDir):
            start=f.rfind("/",0,f.rfind("/")-1)  
        if end==-1:
            sampleName=f[start:]
        else:
            sampleName=f[start:end]
        #strain=f[start:send]
        outSample=outPath+"/"+sampleName
        #if(opts.group):
        #        outSample=outPath+"/"+strain+"_grouped"
        oAligned=outSample+"/aligned.bam"
        oUnmapped=outSample+"/unmapped.bam"
        oSum=outSample+"/summary.txt"
        oSplice=outSample+"/splice_junct.bed"
        if not os.path.exists(outSample):
                os.makedirs(outSample)
        if not os.path.exists(outSample+"/aligned.bam"):
                unmappedR1=""
                unmappedR2=""
                ucount=0
                searchPath=inPath+"/"+sampleName+search
                if(opts.inputSampleDir):
                    searchPath=inPath+"/"+sampleName+search
                #if(opts.group):
                #        searchPath=inPath+"/"+strain+"_*R1*.fq.gz"
                print("search:"+searchPath)
                flist=glob.glob(searchPath)
                for f2 in flist:
                        fr2=f2.replace(opts.pairPrefix+"1",opts.pairPrefix+"2")
                        if(ucount>0):
                                unmappedR1=unmappedR1+","
                                unmappedR2=unmappedR2+","
                        if(opts.inputSampleDir):
                            unmappedR1=f2+"/unmapped.end1.fq.gz"
                            unmappedR2=f2+"/unmapped.end2.fq.gz"
                        else:
                            unmappedR1=unmappedR1+f2
                            unmappedR2=unmappedR2+fr2
                        ucount+=1
                index=indexPath+"/"+opts.indexPref
                if(opts.strainSpecific):
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
                        strain="ACI_EurMcw"
                    elif strain=='Cop':
                        strain="BN"
                    index=indexPath+"/"+strain+opts.indexSuffix
                ##TODO Insert code to substitute a strain name into the indexFile
                arg=str("hisat2 -"+opts.hisatIn+" -p "+str(maxP)+" -k "+str(opts.hisatK))
                if(opts.hisatStrandness):
                        arg=arg+" --rna-strandness "+opts.hisatStrandness
                if(opts.hisatReorder):
                        arg=arg+" --reorder "
                if(opts.hisatSummary):
                        arg=arg+" --new-summary --summary-file "+oSum
                if(opts.hisatSplice):  
                        arg=arg+" --novel-splicesite-outfile "+oSplice
                arg=arg+" -x "+index
                if(opts.paired):
                        arg=arg+" -1 "+unmappedR1+" -2 "+unmappedR2
                else:
                        arg=arg+" -U "+unmappedR1
                arg=arg+" 2> "+outSample+"/err.txt | tee >(samtools view -bS -F 0x4 - > "+oAligned+" 2> /dev/null) >(samtools view -bS -f 0x4 - > "+oUnmapped+" 2> /dev/null) | grep error"
                print("running: "+sampleName)
                print(arg)
                print("\n")
                completed=subprocess.run(arg,shell=True, executable='/bin/bash')
