#!/usr/bin/env python3
import os

from collapse import collapse
from merge import mergeRecord, trimRecord
from trim import trim
from multiprocessing import Process, Pipe
from subprocess import Popen, PIPE
from condense import condense

import BEDtoRef
from sortedcontainers import SortedList
import pysam
import gzip
import sys
import tempfile
import shutil
from configutator import ConfigMap, handleArgs

from valitator import Path, UnsignedInt, PathOrNone, Executable

@ConfigMap(fastq1='fastqs[0]', fastq2='fastqs[1]', bwaArgs='bwa')
def ProDuSe(fastq1: Path, fastq2: Path, reference: str, bed: PathOrNone, output: str, bwaArgs: Executable):
    """

    :param fastq1: Path to the first fastq
    :param fastq2: Path to the second fastq
    :param reference: Path to the reference fasta
    :param bed: Path to a bed file containing regions to restrict reads to (optional)
    :param output:
    :param bcThreshold:
    :param familyThreshold:
    :param sequence:
    :param mask:
    :param bwaArgs:
    :return:
    """
    # Just a heads up, the following code seems out of order. You have to follow the pipes to make sense of it. Good luck Mario ;)
    global config
    temp_dir = tempfile.mkdtemp()

    # Load up bed file coordinates relative to reference index
    if bed:
        coords = BEDtoRef.BEDtoRef(pysam.TabixFile(bed), pysam.TabixFile(reference + ".fai"))
        sys.stderr.write("BED regions loaded.\n")
    else:
        coords = None

    trimOut1Path = os.path.join(temp_dir, 'trimOut1')
    os.mkfifo(trimOut1Path)
    trimOut2Path = os.path.join(temp_dir, 'trimOut2')
    os.mkfifo(trimOut2Path)
    sortInPath = os.path.join(temp_dir, 'sortIn')
    sortOutPath = os.path.join(temp_dir, 'sortOut')
    os.mkfifo(sortInPath)
    os.mkfifo(sortOutPath)
    bwaOutPath = os.path.join(temp_dir, 'bwaOut')
    os.mkfifo(bwaOutPath)

    # Start sort by coord
    sys.stderr.write("Starting sort subprocess..\n")
    sortp = Process(target=pysam.sort, args=("-o", sortOutPath, sortInPath))
    sortp.start()

    #Start clipping and filtering the bwa output
    sys.stderr.write("Starting clipping and read filter subprocess..\n")
    cfp = Process(target=clipAndFilter, args=(bwaOutPath, sortInPath, coords))
    cfp.start()

    # Start bwa subprocess
    if isinstance(bwaArgs, str):
        bwaArgs = bwaArgs.split(' ')
    sys.stderr.write("Starting BWA subprocess..\n")
    bwaArgs += ['-C', reference, trimOut1Path, trimOut2Path]
    bwa = Popen(bwaArgs, stdout=open(bwaOutPath, 'w'))

    # Start trim subprocesses
    sys.stderr.write("Starting trim subprocesses..\n")
    p1Result = []
    p1 = Process(target=callRes, args=(trim, dict(inStream=gzip.open(fastq1), outStream=open(trimOut1Path, 'w'), **config[trim]), p1Result))
    p1.start()
    p2Result = []
    p2 = Process(target=callRes, args=(trim, dict(inStream=gzip.open(fastq2), outStream=open(trimOut2Path, 'w'), **config[trim]), p2Result))
    p2.start()

    sortedBAM = pysam.AlignmentFile(sortOutPath, "r")
    outputBAM = pysam.AlignmentFile(output, "wb", template=sortedBAM)
    sys.stderr.write("Starting collapse subprocess..\n")
    collapseProc = Process(target=collapse, kwargs=dict(inFile=sortedBAM, outFile=outputBAM, **config[collapse]))
    collapseProc.start()

    p1.join()
    p2.join()
    sortp.join()
    cfp.join()
    collapseProc.join()
    bwa.wait()
    shutil.rmtree(temp_dir)

def clipAndFilter(inPath, outPath, regions = None):
    #a = open(inPath, 'r')
    #inPath = '/home/ncm3/bwaout.sam'
    inFile = pysam.AlignmentFile(inPath, 'r', check_sq=False) #TODO: change back to inPath
    outFile = pysam.AlignmentFile(outPath, "w", template=inFile)
    #recordItr = inFile.fetch(until_eof=True)
    while True:
        firstRecord = next(inFile)
        secondRecord = next(inFile)
        # TODO exit loop condition
        if firstRecord.query_name != secondRecord.query_name:
            pass  # TODO bwa was expected to output mate pairs
        if regions and (not BEDtoRef.inCoords(firstRecord.reference_start, regions) \
                or not BEDtoRef.inCoords(firstRecord.reference_end, regions) \
                or not BEDtoRef.inCoords(secondRecord.reference_start, regions) \
                or not BEDtoRef.inCoords(secondRecord.reference_end, regions)):
            continue

        # Condense BWA output
        condense(firstRecord)
        condense(secondRecord)

        forwardRecord = firstRecord if not firstRecord.is_reverse else secondRecord
        reverseRecord = firstRecord if firstRecord.is_reverse else secondRecord

        mergeRecord(reverseRecord, forwardRecord)
        trimRecord(reverseRecord, forwardRecord.reference_end)

        outFile.write(firstRecord)
        outFile.write(secondRecord)

def callRes(func, args, res: list):
    "Process helper to retrieve return value of function"
    res.append(func(**args))

chrList = []
def rangeCmp(a, b):
    '''

    :param a:
    :param b:
    :return:
    '''

    chr1 = a[1]
    if chr1 not in chrList:
        chrList.append(chr1)

    chr2 = b[1]
    if chr2 not in chrList:
        chrList.append(chr2)

    result = chrList.index(chr1) - chrList.index(chr2)

    if result == 0:
        result = a[2] - b[2]

    return result

if __name__ == "__main__":
    ConfigMap(_func='trim')(trim)
    ConfigMap(_func='collapse')(collapse)
    for argmap in handleArgs(sys.argv, (ProDuSe, trim, collapse)):
        global config
        config = argmap
        ProDuSe(**argmap[ProDuSe])
        break