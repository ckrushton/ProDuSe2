import pysam
import sortedcontainers
import math

from CigarIterator import CigarIterator, appendOrInc

class FamilyRecord:
    class Op:
        __slots__ = 'op', 'allele', 'qSum', 'count'
        def __init__(self, op, allele):
            self.op = op
            self.allele = allele
            self.maxQ = 0
            self.qSum = 0
            self.count = 0

        def __iadd__(self, other: int):
            if self.maxQ < other:
                self.maxQ = other
            self.qSum += other
            self.count += 1

        def merge(self, op):
            if self.maxQ < op.maxQ:
                self.maxQ = op.maxQ
            self.qSum += op.qSum
            self.count += op.count

        def inSeq(self) -> bool:  # Returns true if the operation has a sequence coordinate
            return self.op in (pysam.CMATCH, pysam.CINS, pysam.CSOFT_CLIP, pysam.CEQUAL, pysam.CDIFF)

    __slots__ = 'name', 'pos', 'maxMapQ', 'cols', 'size'
    def __init__(self, name: str, pos: int, record: pysam.AlignedSegment = None):
        self.name = name
        self.pos = pos
        self.cols = []
        self.maxMapQ = 0
        self.size = 0
        if record:
            self.aggregate(record)

    def getOpsAt(self, pos: int):
        def OPitr():
            for op in self.cols[pos]:
                yield op

        if len(self.cols) > pos:
            return OPitr
        else:
            return iter(())

    # Assumes first unclipped op in recordItr aligns to opCount[0]
    def aggregate(self, record: pysam.AlignedSegment):
        if self.maxMapQ < record.mapping_quality:
            self.maxMapQ = record.mapping_quality

        recordItr = CigarIterator(record)
        recordItr.skipClipped()
        i = 0
        while recordItr.next():
            if len(self.cols) < i:
                op = self.Op(recordItr.op, recordItr.seqBase)
                op += recordItr.baseQual or 0
                self.cols.append(sortedcontainers.SortedListWithKey(key=lambda x: (x.op, x.allele), iterable=[op]))
            else:
                self.cols[i][(recordItr.op, recordItr.seqBase)] += recordItr.baseQual or 0
            i += 1
        self.size += 1

    def __iadd__(self, other):
        if self.maxMapQ < other.maxMapQ:
            self.maxMapQ = other.maxMapQ
        pos = 0
        # Merge in other.cols
        while pos < len(self.cols) and pos < len(other.ops):
            for k, op in other.ops[pos]:
                self.cols[pos][k].merge(op)
            pos += 1

        # If other.cols is longer, extend self.cols
        while pos < len(other.ops):
            self.cols.append(other.ops[pos])
            pos += 1

    def __len__(self):
        return self.size

    @property
    def is_positive_strand(self) -> bool:

        pass #TODO

    def toPysam(self) -> pysam.AlignedSegment:
        seq = ""
        qual = []
        ops = []
        fQ = []
        fC = []
        for col in self.cols:
            bestOp = col[0]
            totalN = 0
            totalQ = 0
            for op in col:
                totalN += op.count
                totalQ += op.qSum
                if bestOp.qSum < op.qSum:
                    bestOp = op
            if bestOp.inSeq():
                seq += bestOp.allele
                qual += self.mapQ
                appendOrInc(ops, [op.op, 1])
                fQ += [int(round(-10*math.log10((totalQ - bestOp.qSum) * bestOp.count / (totalQ * totalN))))]
                fC += [totalN]
        record = pysam.AlignedSegment()
        record.query_name = self.name
        record.reference_start = self.pos
        record.mapping_quality = self.maxMapQ
        record.cigartuples = ops
        record.query_sequence = seq
        record.query_qualities = qual
        # Added tags
        # fQ:B:C Integer array containing Phred score of wrong base chosen during collapse
        # fC:B:I Integer array containing pairs (simmilar to a CIGAR) of count and depth representing, from the start of the alignment, the depth of the family at a position
        record.set_tags([('fQ', fQ), ('fC', fC)])
        return record