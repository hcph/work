#!/usr/bin/env python
#coding: utf8

###### Import Modules
import sys, os
import re
import glob
import gzip
# import json
# import pprint
# import cPickle as pickle

###### Document Decription
'''  '''

###### Version and Date
prog_version = '1.1.0'
prog_date = '2015-04-02'

###### Usage
usage = '''

     Version %s  by Vincent-Li  %s

     Usage: %s <barcodeList> <fastqFile> <outDir> 
''' % (prog_version, prog_date, os.path.basename(sys.argv[0]))

######## Global Variable


#######################################################################
############################  BEGIN Class  ############################
#######################################################################


##########################################################################
############################  BEGIN Function  ############################
##########################################################################
def createBarcodeDict(fileName, barh, ofhDict, outd, err=1):
    barSet = set()
    ########### Open and read input file
    try:
        fhIn = open(fileName, 'r')
    except Exception as e:
        raise e
    
    for line in fhIn:
        info = line.split()
        if len(info) < 2:
            continue
        if re.match(r'^\d+$', info[0]):
            info[0] = "barcode_" + info[0]
        createBarHash(info[0], info[1], barh, err)

        barSet.add(info[1])

        ########### Open and write file
        try:
            newFh = open(os.path.join(outd, "%s.fq" % info[0]), 'w')
        except Exception, e:
            raise e
        ofhDict[info[0]] = newFh;

        ##### ambiguous.fq
        if "ambiguous" in barh.values():
            try:
                ambiFh = open(os.path.join(outd, "%s.fq" % "ambiguous"), 'w')
            except Exception, e:
                raise e
            ofhDict["ambiguous"] = ambiFh;

    fhIn.close
    return barSet

def createBarHash(name, seq, bh, err=1):
    if err == 0:
        bh[seq] = name
        return bh
    baseList = ('A', 'C', 'G', 'T', 'N')
    for idx in xrange(0, len(seq)):
        for b in baseList:
            tmpSeq = seq[:idx] + b + seq[idx+1:]
            if err > 1:
                createBarHash(name, tmpSeq, bh, err-1)
           
            if tmpSeq in bh and bh[tmpSeq] != name:
                bh[tmpSeq] = 'ambiguous'
            else:
                bh[tmpSeq] = name
    return bh

def splitBarcode(fqFile, barh, fhd, fc=36, lc=45):
    statDict = {}
    ########### Open and read input file
    try:
        if fqFile.endswith('.gz'):
            fqIn = gzip.open(fqFile, 'r')
        else:
            fqIn = open(fqFile, 'r')
    except Exception as e:
        raise e
    
    idLine = "1"
    while idLine:
        idLine = fqIn.readline()
        seqLine = fqIn.readline()
        plusLine = fqIn.readline()
        qualLine = fqIn.readline()
        barcodeSeq = seqLine[fc-1:lc]

        if barcodeSeq not in statDict:
            statDict[barcodeSeq] = 0
        statDict[barcodeSeq] += 1

        if barcodeSeq in barh:
            ofh = fhd[barh[barcodeSeq]]
            #print ofh
            ofh.write(idLine)
            ofh.write(seqLine[:fc-1] + seqLine[lc:])
            ofh.write(plusLine)
            ofh.write(qualLine[:fc-1] + qualLine[lc:])

    fqIn.close
    return statDict

def outStat(statd, barh, bset, outDir):
    total = sum(statd.values())
    klist = sorted(statd.keys(), key=lambda x: statd[x], reverse=True)
    statDict = {}
    ########### Open and write file
    try:
        tagStat = open(os.path.join(outDir, "tagStatistics.txt"), 'w')
        splitStat = open(os.path.join(outDir, "splitRate.txt"), 'w')
    except Exception, e:
        raise e
    
    for k in klist:
        bar = "unknown"
        if k in barh:
            bar = barh[k]
            if bar not in statDict:
                statDict[bar] = {"correct": 0, "corrected": 0}
            if k in bset:
                statDict[bar]['correct'] += statd[k]
            else:
                statDict[bar]['corrected'] += statd[k]
        tagStat.write("%s\t%s\t%s\t%.3f\n" % (k ,bar, statd[k], 100.0*statd[k]/total))
    tagStat.close
    splitStat.write("#barcode\tCorrect\tCorrected\tTotal\tPct%\n")
    allCorrect = 0
    allCorrected = 0
    for b in sorted(statDict.keys()):
        splitedNum = statDict[b]['correct'] + statDict[b]['corrected']
        splitStat.write("%s\t%s\t%s\t%s\t%.3f\n" % (b, statDict[b]['correct'], statDict[b]['corrected'], splitedNum, 100.0*splitedNum/total))
        allCorrect += statDict[b]['correct']
        allCorrected += statDict[b]['corrected']
    splitStat.write("Total\t%s\t%s\t%s\t%.3f\n" % (allCorrect, allCorrected, allCorrect + allCorrected, 100.0*(allCorrect + allCorrected)/total))

######################################################################
############################  BEGIN Main  ############################
######################################################################
#################################
##
##   Main function of program.
##
#################################
def main():
    
    ######################### Phrase parameters #########################
    import argparse
    ArgParser = argparse.ArgumentParser(usage = usage, version = prog_version)
    ArgParser.add_argument("-e", "--errNum", action="store", dest="errNum", type=int, default=2, metavar="INT", help=" Allow mismatch count in barcode sequence. [%(default)s]")
    ArgParser.add_argument("-f", "--firstCycle", action="store", dest="firstCycle", type=int, default=36, metavar="INT", help="First cylce of barcode. [%(default)s]")
    ArgParser.add_argument("-l", "--lastCycle", action="store", dest="lastCycle", type=int, default=45, metavar="INT", help="Last cycle of barcode. [%(default)s]")

    (para, args) = ArgParser.parse_known_args()

    if len(args) != 3:
        ArgParser.print_help()
        print >>sys.stderr, "\nERROR: The parameters number is not correct!"
        sys.exit(1)
    else:
        (barcodeList, fastqFile, outDir) = args

    ############################# Main Body #############################
    barcodeHash = {}
    outFh = {}
    barcodeSet = createBarcodeDict(barcodeList, barcodeHash, outFh, outDir, para.errNum)
    splitStat = splitBarcode(fastqFile, barcodeHash, outFh, para.firstCycle, para.lastCycle)
    outStat(splitStat, barcodeHash, barcodeSet, outDir)


#################################
##

##   Start the main program.
##
#################################
if __name__ == '__main__':
    main()

################## God's in his heaven, All's right with the world. ##################
