import pysam
import math

try:
    from CigarIterator import CigarIterator, appendOrInc
# If installed
except ModuleNotFoundError:
    from ProDuSe.CigarIterator import CigarIterator, appendOrInc

class FamilyRecord:
    class Op:
        __slots__ = 'op', 'allele', 'qSum', 'maxQ', 'count'
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
            return self

        def merge(self, op):
            if self.maxQ < op.maxQ:
                self.maxQ = op.maxQ
            self.qSum += op.qSum
            self.count += op.count

        def inSeq(self) -> bool:  # Returns true if the operation has a sequence coordinate
            return self.op in (pysam.CMATCH, pysam.CINS, pysam.CSOFT_CLIP, pysam.CEQUAL, pysam.CDIFF)

    __slots__ = 'name', 'pos', 'maxMapQ', 'cols', 'members'
    def __init__(self, name: str, pos: int, record: pysam.AlignedSegment = None):
        self.name = name
        self.pos = pos
        self.cols = []
        self.maxMapQ = 0
        self.members = []
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
        try:
            startPos = record.get_tag('OS')  # type: int
        except KeyError:
            startPos = record.reference_start
        recordItr = CigarIterator(record)
        recordItr.skipClipped()
        i = startPos - record.reference_start
        while recordItr.valid:
            if len(self.cols) <= i:
                op = self.Op(recordItr.op, recordItr.seqBase)
                op += recordItr.baseQual or 0
                pos = {}
                pos[(op.op, op.allele)] = op
                self.cols.append(pos)
            else:
                op = self.cols[i].get((recordItr.op, recordItr.seqBase))
                if op:
                    op += recordItr.baseQual or 0
                else:
                    self.cols[i][(recordItr.op, recordItr.seqBase)] = self.Op(recordItr.op, recordItr.seqBase)
            i += 1
            recordItr.next()
        self.members.append(record.query_name)

    def __iadd__(self, other: 'FamilyRecord'):
        if self.maxMapQ < other.maxMapQ:
            self.maxMapQ = other.maxMapQ
        pos = 0
        # Merge in other.cols
        while pos < len(self.cols) and pos < len(other.cols):
            for k, op in other.cols[pos].items():
                if k in self.cols[pos]:
                    self.cols[pos][k].merge(op)
                else:
                    self.cols[pos][k] = op
            pos += 1

        # If other.cols is longer, extend self.cols
        while pos < len(other.cols):
            self.cols.append(other.cols[pos])
            pos += 1
        self.members += other.members
        return self

    def __len__(self):
        return len(self.members)

    @property
    def is_positive_strand(self) -> bool:

        return False #TODO

    def toPysam(self) -> pysam.AlignedSegment:
        seq = ""
        qual = []
        ops = []
        fQ = []
        fC = []

        # Compile sequence and added tags
        # fQ:B:C Integer array containing Phred score of wrong base chosen during collapse
        # fC:B:I Integer array containing pairs (simmilar to a CIGAR) of count and depth representing, from the start of the alignment, the depth of the family at a position
        for col in self.cols: #dict
            colOps = list(col.values())
            bestOp = colOps[0]
            totalN = 0
            totalQ = 0
            for op in colOps:
                totalN += op.count
                totalQ += op.qSum
                if bestOp.qSum < op.qSum:
                    bestOp = op
            if bestOp.inSeq():
                seq += bestOp.allele
                qual.append(bestOp.maxQ)
                appendOrInc(ops, [op.op, 1])
                if (totalQ - bestOp.qSum) > 0:
                    fQ.append(int(round(-10*math.log10((totalQ - bestOp.qSum) * bestOp.count / (totalQ * totalN)))))
                else:
                    fQ.append(93) # Max Phred Score
                fC.append(totalN)
        tags = []
        if fQ:
            tags.append(('fQ', fQ))
        if fC:
            tags.append(('fC', fC))

        # Copy into pysam record
        record = pysam.AlignedSegment()
        record.query_name = self.pos + self.name
        record.reference_start = self.pos
        record.mapping_quality = self.maxMapQ
        record.flag = 17
        record.cigartuples = ops
        record.query_sequence = seq
        record.query_qualities = qual
        if tags:
            record.set_tags(tags)
        return record