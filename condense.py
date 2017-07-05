import pysam
from CigarIterator import CigarIterator, appendOrInc

def condense(record: pysam.AlignedSegment):

    if record.reference_length == 0:
        # No work needs to be done
        return
    seq = ""
    qual = []
    ops = []
    itr = CigarIterator(record)
    clipped = itr.skipClipped()
    if clipped:
        ops.append((pysam.CHARD_CLIP, clipped))
    lastPos = itr.opPos
    while itr.skipToNonRef(): # TODO: check for MD tag and add if not present
        if itr.inSeq:
            seq += itr.seqBase
            qual.append(itr.baseQual)
        dist = itr.opPos - lastPos
        if dist > 0:
            appendOrInc(ops, [pysam.CREF_SKIP, dist])
        if itr.op == pysam.CMATCH:
                appendOrInc(ops, [pysam.CDIFF, 1])
        else:
            appendOrInc(ops, [itr.op, 1])
    record.query_sequence = seq
    record.query_qualities = qual
    record.cigartuples = ops

if __name__ == "__main__":
    pass #TODO