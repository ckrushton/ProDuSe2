import pysam
from sys import maxsize

from CigarIterator import CigarIterator, appendOrInc

default_cost = {
    pysam.CMATCH: lambda x: -x,
    pysam.CEQUAL: lambda x: -x,
    pysam.CDIFF: lambda x: -x,
    pysam.CREF_SKIP: lambda x: -x,
    pysam.CINS: lambda x: 6 + 1*(x-1),
    pysam.CDEL: lambda x: 3 + 1*(x-1)
}

def trimRecord(record: pysam.AlignedSegment, start: int = 0, end: int = maxsize):
    ops = []
    for i in CigarIterator(record):
        if (i.refPos < start or i.refPos > end) and not i.clipped():
            if i.inSeq():
                appendOrInc(ops, [pysam.CSOFT_CLIP, 1])
            else:
                appendOrInc(ops, [pysam.CHARD_CLIP, 1])
        else:
            appendOrInc(ops, [i.getOp(), 1])
    record.cigartuples = ops

def mergeRecord(fromRecord: pysam.AlignedSegment, toRecord: pysam.AlignedSegment, refStart: int = -1, refEnd: int = maxsize, costs = default_cost) -> list:
    stats = []
    seq = ""
    qual = []
    ops = []
    toItr = CigarIterator(toRecord)
    fromItr = CigarIterator(fromRecord)

    rStart = fromRecord.reference_start if fromRecord.reference_start < toRecord.referece_start else toRecord.reference_start
    rEnd = (fromRecord.reference_end if fromRecord.reference_end < toRecord.referece_end else toRecord.reference_end) - 1
    if refStart < rStart:
        refStart = rStart
    if refEnd > rEnd:
        refEnd = rEnd

    toOptimal = None

    while toItr.refPos <= refEnd or fromItr.refPos <= refEnd:
        if toItr.getOp() == fromItr.getOp():
            if toItr.inSeq():  #.getOp() in [pysam.CMATCH, pysam.CEQUAL, pysam.CDIFF]:
                if toItr.getBaseQual() > fromItr.getBaseQual():
                    seq += toItr.getSeqBase()
                    qual += [toItr.getBaseQual()]
                else:
                    seq += fromItr.getSeqBase()
                    qual += [fromItr.getBaseQual()]
            appendOrInc(ops, [toItr.getOp(), 1])
            continue
        if toOptimal == None:
            # Dont calculate costs if unnecessary 
            toCost = 0
            fromCost = 0
            tItr = CigarIterator(toRecord)
            fItr = CigarIterator(fromRecord)
            tItr.skipToRefPos(refStart)
            fItr.skipToRefPos(refStart)
            while True:
                toCost += costs[tItr.getOp()](tItr.opEnd() - tItr.opStart + 1)
                if not tItr.nextOp() or tItr.clipped(): break

            while True:
                fromCost += costs[fItr.getOp()](fItr.opEnd() - fItr.opStart + 1)
                if not fItr.nextOp() or fItr.clipped(): break

            toOptimal = toCost > fromCost
        elif toOptimal:
            if toItr.inSeq():
                seq += toItr.getSeqBase()
                qual += [toItr.getBaseQual()]
            appendOrInc(ops, [toItr.getOp(), 1])
        else:
            if fromItr.inSeq():
                seq += fromItr.getSeqBase()
                qual += [fromItr.getBaseQual()]
            appendOrInc(ops, [fromItr.getOp(), 1])

    toRecord.cigartuples = ops
    toRecord.query_sequence = seq
    toRecord.query_qualities = qual
