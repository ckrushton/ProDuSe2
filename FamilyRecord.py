import pysam
import math

from CigarIterator import CigarIterator, appendOrInc

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

    __slots__ = 'name', 'pos', 'maxMapQ', 'cols', 'members', 'mate', 'barcode', 'refID', 'forwardCounter'
    def __init__(self, name: str, pos: int, barcode: str, record: pysam.AlignedSegment):
        self.name = name
        self.pos = pos
        self.cols = []
        self.maxMapQ = 0
        self.members = []
        self.mate = None #type: FamilyRecord
        self.barcode = barcode
        self.refID = record.reference_id
        self.forwardCounter = 0
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
        self.forwardCounter += 1 - (2 * record.is_reverse) # Increment if forward, decrement if reverse

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
        self.forwardCounter += other.forwardCounter
        return self

    def __len__(self):
        return len(self.members)

    def getConcensusAt(self, pos: int) -> (int, str or None, int or None, int or None, int or None):
        """
        Compute the consensus operation at a specified position in the record.
        Everything but the returned op may be None if not relevant to op.
        :param pos: The position in the record in cigar space
        :return: Tuple(op code, base, qScore, depth at pos, Phred score of wrong op chosen)
        """
        colOps = list(self.cols[pos].values())
        bestOp = colOps[0]
        totalN = 0
        totalQ = 0
        base = None
        qual = None
        fC = None
        fQ = None
        for op in colOps:
            totalN += op.count
            totalQ += op.qSum
            if bestOp.qSum < op.qSum:
                bestOp = op
        if bestOp.inSeq():
            base = bestOp.allele
            qual = bestOp.maxQ
            op = op.op
            if (totalQ - bestOp.qSum) > 0:
                fQ = int(round(-10 * math.log10((totalQ - bestOp.qSum) * bestOp.count / (totalQ * totalN))))
            else:
                fQ = 93  # Max Phred Score
            fC = totalN
        return op, base, qual, fC, fQ

    def __iter__(self):
        """
        Iterate every position returning result of self.getConcensusAt()
        :return: Generator that iterates all positions in self.cols
        """
        for pos in range(len(self.cols)):
            yield self.getConcensusAt(pos)