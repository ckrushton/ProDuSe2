from collapse import collapse
from merge import mergeRecord, trimRecord
from trim import trim
from multiprocessing import Process, Pipe
from subprocess import Popen, PIPE
from condense import condense
import pysam

def callRes(func, args, res):
    "Process helper to retrieve return value of function"
    res = func(**args)

def __init__():

    # Load args

    trimIn1, trimOut1 = Pipe(False)
    trimIn2, trimOut2 = Pipe(False)

    # Start trim subprocesses
    p1Result = []
    p1 = Process(target=callRes, args=(trim, (open(fastq1), open(trimOut1, 'w'), threshold, sequence, reverse), p1Result))
    p1.start()
    p2Result = []
    p2 = Process(target=callRes, args=(trim, (open(fastq2), open(trimOut2, 'w'), threshold, sequence, reverse), p2Result))
    p2.start()

    # Start bwa subprocess
    bwaArgs = ['bwa', 'mem', '-C', reference, '<&' + str(trimIn1.fileno()), '<&' + str(trimIn2.fileno())]
    bwa = Popen(bwaArgs, stdout=PIPE, pass_fds=(trimIn1.fileno(), trimIn2.fileno()))

    #Start sort by coord
    sortIn, sortOut = Pipe(False)
    pysam.sort("-o", '<&' + sortOut.fileno(), '<&' + sortIn.fileno())

    unsortedBAM = pysam.AlignmentFile(sortIn, "wb")
    sortedBAM = pysam.AlignmentFile(sortOut, "rb")

    outputBAM = pysam.AlignmentFile(outputPath, "wb")

    collapseProc = Process(target=collapse, args=(sortedBAM, outputBAM, threshold, mask))
    collapseProc.start()

    recordItr = pysam.AlignmentFile(bwa.stdout).fetch(until_eof=True)
    while True:
        firstRecord = recordItr.next()
        secondRecord = recordItr.next()
        # TODO exit loop condition
        if firstRecord.query_name != secondRecord.query_name:
            pass # TODO bwa was expected to output mate pairs

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