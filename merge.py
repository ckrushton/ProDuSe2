import pysam
from sys import maxsize
from CigarIterator import CigarIterator

default_cost = {
    pysam.CMATCH: lambda x: -x,
    pysam.CREF_SKIP: lambda x: -x,
    pysam.CINS: lambda x: 6 + 1*(x-1),
    pysam.CDEL: lambda x: 3 + 1*(x-1)
}

def appendOrInc(ops: [], op: list):
    if len(ops) > 0 and ops[-1][0] == op[0]:
        ops[-1][1] += op[1]
    else:
        ops.append(op)

def mergeRecord(fromRecord: pysam.AlignedSegment, toRecord: pysam.AlignedSegment, refStart: int = -1, refEnd: int = maxsize, costs = default_cost) -> list:
    stats = []
    seq = toRecord.query_sequence
    qual = toRecord.query_qualities
    ops = []
    toItr = CigarIterator(toRecord)
    fromItr = CigarIterator(fromRecord)

    rStart = fromRecord.reference_start if fromRecord.reference_start < toRecord.referece_start else toRecord.reference_start
    rEnd = (fromRecord.reference_end if fromRecord.reference_end < toRecord.referece_end else toRecord.reference_end) - 1
    if refStart < rStart:
        refStart = rStart
    if refEnd > rEnd:
        refEnd = rEnd

    toItr.skipToRefPos(refStart)
    fromItr.skipToRefPos(refStart)

    toCost = 0
    fromCost = 0

    while True:
        toCost = costs[toItr.getOp()](toItr.opEnd() - toItr.opStart + 1)
        if not toItr.nextOp() or toItr.clipped(): break
    
    while True:
        fromCost = costs[fromItr.getOp()](fromItr.opEnd() - fromItr.opStart + 1)
        if not fromItr.nextOp() or fromItr.clipped(): break

    toOptimal = toCost > fromCost

    while toItr.refPos <= refEnd or fromItr.refPos <= refEnd:
        if toItr.getOp() == fromItr.getOp():
            if toItr.getOp() in [pysam.CMATCH, pysam.CEQUAL, pysam.CDIFF]:

