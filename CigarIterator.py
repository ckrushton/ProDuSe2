import pysam
import re

MDOps = r'(\d*)\^?([A-Za-z])'

def appendOrInc(ops: [], op: list):
    if len(ops) > 0 and ops[-1][0] == op[0]:
        ops[-1][1] += op[1]
    else:
        ops.append(op)

class CigarIterator(object):
    def __init__(self, record: pysam.AlignedSegment):
        self.record = record
        self.ops = record.cigartuples  # List of CIGAR operations
        self.md = []  # Reference bases from MD tag
        self._buildMD()

    def rewind(self):
        self.opsI = 0  # Current index in ops
        self.opPos = -1  # Current position in operation including all previous operations
        self.opStart = 0  # Start of current operation
        self.seqPos = 0  # Current sequence position
        self.refStart = self.record.reference_start  # Aligned starting position of unclipped sequence
        self.refPos = self.refStart  # Current reference position
        self.mdI = 0  # Current index in MD

    def __iter__(self):
        self.rewind()
        return self

    def __next__(self):
        self.opPos += 1
        if self.opEnd() <= self.opPos:
            self.opStart += self.opLength()
            self.opsI += 1
        if self.opsI >= len(self.ops):
            raise StopIteration
        while self._getMD() and self._getMD()[1] < self.opPos:
            self.mdI += 1
        if self.inSeq():
            self.seqPos += 1
        if self.inRef():
            self.refPos += 1

    def next(self) -> bool:
        try:
            self.__next__()
            return True
        except StopIteration:
            return False

    def nextOp(self) -> bool:
        self.opsI += 1
        if self.opsI >= len(self.ops):
            return False
        dist = self.opEnd() - self.opPos
        self.opStart += self.opLength()
        self.opPos = self.opStart
        while self._getMD() and self._getMD()[1] < self.opPos:
            self.mdI += 1
        if self.inSeq():
            self.seqPos += dist
        if self.inRef():
            self.refPos += dist
        return True

    def _buildMD(self):
        if self.record.has_tag("MD"):
            mdStr = self.record.get_tag("MD")
            mdOffset = self.record.query_alignment_start # Cigar pos will correlate to number of clipped bases
            for mdOp in re.finditer(MDOps, mdStr):
                matchCount, refBase = mdOp.group(1, 2)
                if matchCount != "":
                    mdOffset += int(matchCount)
                self.md.append((refBase, mdOffset))
                mdOffset += 1

    def _getMD(self) -> tuple:
        return self.md[self.mdI] if self.mdI >= 0 and self.mdI < len(self.md) else None

    def _getOpMD(self):
        md = self._getMD()
        return md if md and md[1] == self.opPos else None

    def opLength(self) -> int:
        return self.ops[self.opsI][1] if len(self.ops) else 0

    def opEnd(self) -> int:
        l = self.opLength()
        if l == 0:
            return self.opStart
        else:
            return self.opStart + l - 1

    def skipClipped(self) -> int:
        count = 0
        if self.opLength() > 0:
            while self.opsI < len(self.ops) and self.clipped():
                count += self.opLength()
                self.opsI += 1
            self.opStart = count
            self.opPos = self.opStart

        return count

    def skipToRefPos(self, pos: int):
        pass #TODO

    def skipToNonRef(self) -> bool: # Move iterator to next non-reference cigar position (variant in MD tag)
        if len(self.md) == 0:
            return False
        if self.md[self.mdI][1] == self.opPos:
            self.mdI += 1
        if self.mdI >= len(self.md) or self.opsI >= len(self.ops):
            return False
        newPos = self.md[self.mdI][1]
        if self.opPos < 0:
            self.opPos = 0
        dist = newPos - self.opPos

        while dist > 0:
            if newPos < self.opEnd():
                if self.inSeq():
                    self.seqPos += dist
                if self.inRef():
                    self.refPos += dist
                dist = 0
                self.opPos = newPos
            else:
                d = self.opEnd() - self.opPos
                dist -= d
                if self.inSeq():
                    self.seqPos += d
                if self.inRef():
                    self.refPos += d
                self.opStart += self.opLength()
                self.opPos = self.opStart
                self.opsI += 1
                if self.opsI >= len(self.ops):
                    return False
        return True



    def inRef(self) -> bool: # Returns true if the passed operation has a reference coordinate
        return self.getOp() in (pysam.CMATCH, pysam.CDEL, pysam.CREF_SKIP, pysam.CEQUAL, pysam.CDIFF)

    def inSeq(self) -> bool: # Returns true if the passed operation has a sequence coordinate
        return self.getOp() in (pysam.CMATCH, pysam.CINS, pysam.CSOFT_CLIP, pysam.CEQUAL, pysam.CDIFF)

    def clipped(self):
        return self.getOp() in (pysam.CHARD_CLIP, pysam.CSOFT_CLIP)

    def getRefBase(self) -> str:
        return (self.getSeqBase() if self.matchesRef() else self.md[self.mdI][0]) if self.inRef() else ""

    def matchesRef(self) -> bool:
        return len(self.md) == 0 or self.md[self.mdI][1] != self.opPos #self.getSeqBase() == self.getRefBase()

    def getSeqBase(self) -> str:
        return self.record.query_alignment_sequence[self.seqPos] if self.inSeq() else ""

    def setSeqBase(self, str) -> bool:
        if self.inSeq():
            self.record.query_alignment_sequence[self.seqPos] = str
        else:
            return False
        return True

    def getBaseQual(self) -> int:
        return self.record.query_alignment_qualities[self.seqPos] if self.inSeq() else None

    def setBaseQual(self, qual: int) -> bool:
        if self.inSeq():
            self.record.query_alignment_qualities[self.seqPos] = qual
        else:
            return False
        return True

    def getOp(self) -> int:
        return self.ops[self.opsI][0]