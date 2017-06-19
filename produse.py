import collapse
from merge import mergeRecord, trimRecord
from trim import trim
import trim
import filter
from multiprocessing import Process, Pipe
from subprocess import Popen, PIPE
from condense import condense
import BEDtoRef
from sortedcontainers import SortedList
import pysam
import os
from configutator import ValidateArgs, Assert

@ValidateArgs(os.path.isfile, os.path.isfile, os.path.isfile, Assert(lambda x: not x or os.path.isfile(x), FileNotFoundError("The bed file passed to {func.__name__} does not exist: {val}")), lambda x: x >= 0, lambda x: x >= 0, lambda x,args: len(x) == len(args["mask"]))
def ProDuSe(fastq1, fastq2, reference, bed, output, bcThreshold, familyThreshold, sequence, mask, bwaArgs):
    # Just a heads up, the following code seems out of order. You have to follow the pipes to make sense of it. Good luck Mario ;)

    # Load up bed file coordinates relative to reference index
    if bed:
        coords = BEDtoRef.BEDtoRef(pysam.TabixFile(bed), pysam.TabixFile(reference + ".fai"))
    else:
        coords = SortedList()

    trimIn1, trimOut1 = Pipe(False)
    trimIn2, trimOut2 = Pipe(False)

    # Start trim subprocesses
    p1Result = []
    p1 = Process(target=callRes, args=(trim, (open(fastq1), open(trimOut1, 'w'), bcThreshold, sequence), p1Result))
    p1.start()
    p2Result = []
    p2 = Process(target=callRes, args=(trim, (open(fastq2), open(trimOut2, 'w'), bcThreshold, sequence), p2Result))
    p2.start()

    # Start bwa subprocess
    bwaArgs += ['<&' + str(trimIn1.fileno()), '<&' + str(trimIn2.fileno())]
    bwa = Popen(bwaArgs, stdout=PIPE, pass_fds=(trimIn1.fileno(), trimIn2.fileno()))

    #Start sort by coord
    sortIn, sortOut = Pipe(False)
    pysam.sort("-o", '<&' + sortOut.fileno(), '<&' + sortIn.fileno())

    unsortedBAM = pysam.AlignmentFile(sortIn, "wb")
    sortedBAM = pysam.AlignmentFile(sortOut, "rb")

    outputBAM = pysam.AlignmentFile(output, "wb")

    collapseProc = Process(target=collapse.collapse, args=(sortedBAM, outputBAM, familyThreshold, mask))
    collapseProc.start()

    recordItr = pysam.AlignmentFile(bwa.stdout).fetch(until_eof=True)
    while True:
        firstRecord = recordItr.next()
        secondRecord = recordItr.next()
        # TODO exit loop condition
        if firstRecord.query_name != secondRecord.query_name:
            pass # TODO bwa was expected to output mate pairs
        if bed and not BEDtoRef.inCoords(firstRecord.reference_start, coords) or not BEDtoRef.inCoords(firstRecord.reference_start, coords):
            continue

        # Condense BWA output
        condense(firstRecord)
        condense(secondRecord)

        forwardRecord = firstRecord if not firstRecord.is_reverse else secondRecord
        reverseRecord = firstRecord if firstRecord.is_reverse else secondRecord

        mergeRecord(reverseRecord, forwardRecord)
        trimRecord(reverseRecord, forwardRecord.reference_end)

        unsortedBAM.write(firstRecord)
        unsortedBAM.write(secondRecord)

    p1.join()
    p2.join()
    collapseProc.join()

def callRes(func, args, res):
    "Process helper to retrieve return value of function"
    res = func(**args)

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

def getArgs(parser):
    # Universal Args
    pass

if __name__ == "__main__":
    # Look, we both know this is terrible. But resolve_conflicts is also terrible. Thus, I don't have any choice but to copy-paste arguments over

    bwaArgs = ['bwa', 'mem', '-C', reference]


    getArgs(parser)
    trim.getArgs(parser)
    collapse.getArgs(parser)
    filter.getArgs(parser)
