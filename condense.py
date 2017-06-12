import pysam
from CigarIterator import CigarIterator

def appendOrInc(ops: [], op: list):
    if len(ops) > 0 and ops[-1][0] == op[0]:
        ops[-1][1] += op[1]
    else:
        ops.append(op)

def condense(inFile: pysam.AlignmentFile, outFile: pysam.AlignmentFile):
    for record in inFile.fetch(until_eof=True):
        seq = ""
        qual = []
        ops = []
        itr = CigarIterator(record)
        clipped = itr.skipClipped()
        if clipped:
            ops.append((pysam.CHARD_CLIP, clipped))
        lastPos = itr.opPos
        while itr.skipToNonRef(): # TODO: check for MD tag and add if not present
            if itr.inSeq():
                seq += itr.getSeqBase()
                qual.append(itr.getBaseQual())
            dist = itr.opPos - lastPos
            if dist > 0:
                appendOrInc(ops, [pysam.CREF_SKIP, dist])
            if itr.getOp() == pysam.CMATCH:
                    appendOrInc(ops, [pysam.CDIFF, 1])
            else:
                appendOrInc(ops, [itr.getOp(), 1])
        record.query_sequence = seq
        record.query_qualities = qual
        record.cigartuples = ops
        outFile.write(record)