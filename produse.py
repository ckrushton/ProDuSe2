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

    # Condense BWA output into BAM
    condense(pysam.AlignmentFile(bwa.stdout), pysam.AlignmentFile(sortIn))

    

    p1.join()
    p2.join()