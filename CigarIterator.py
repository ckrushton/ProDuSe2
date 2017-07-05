import pysam
import re

MDOps = r'(\d*)\^?([A-Za-z])'

def appendOrInc(ops: list, op: list):
    if op[1] <= 0:
        return
    if len(ops) > 0 and ops[-1][0] == op[0]:
        ops[-1][1] += op[1]
    else:
        ops.append(list(op))

class CigarIterator(object):
    def __init__(self, record: pysam.AlignedSegment):
        self.record = record
        self.ops = record.cigartuples or []  # List of CIGAR operations
        self.md = None  # Reference bases from MD tag
        self.rewind()

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
        if self.next():
            return self
        else:
            raise StopIteration

    @property
    def valid(self):
        return self.opsI < len(self.ops)

    def step(self, i: int):
        #TODO support negative step
        if i < 0: raise NotImplementedError("Negative stepping not yet supported.")
        if self.opPos < 0:
            self.opPos = 0
        self.opPos += i
        if not self.valid:
            return False
        while self.opEnd < self.opPos:
            delta = self.opEnd - self.opPos + 1
            i -= delta
            if self.inSeq:
                self.seqPos += delta
            if self.inRef:
                self.refPos += delta
            self.opStart += self.opLength
            self.opsI += 1
            if not self.valid:
                return False
        if self.inSeq:
            self.seqPos += i
        if self.inRef:
            self.refPos += i

        return True

    def next(self) -> bool:
        if self.opPos < 0:
            self.opPos = 0
            return self.valid
        return self.step(1)

    #def prev(self) -> bool:
    #    return self.step(-1)

    def nextOp(self) -> bool:
        if self.opPos < 0:
            self.opPos = 0
            return self.valid
        if not self.valid or not len(self.ops):
            return False
        dist = self.ops[self.opsI][1] - (self.opPos - self.opStart)
        self.opStart += self.opLength
        self.opPos = self.opStart
        if self.inSeq:
            self.seqPos += dist
        if self.inRef:
            self.refPos += dist
        self.opsI += 1
        return self.valid

    def _buildMD(self):
        if self.record.has_tag("MD"):
            self.md = []
            mdStr = self.record.get_tag("MD")
            i = 0
            pos = 0
            mdOffset = self.record.query_alignment_start # Cigar pos will correlate to number of clipped bases
            for mdOp in re.finditer(MDOps, mdStr):
                matchCount, refBase = mdOp.group(1, 2)
                mdOffset += int(matchCount or 0)
                while pos <= mdOffset: # Scan CIGAR for insertions and add to offset as MD does not include insertions in MD coordinate space
                    pos += self.ops[i][1]
                    if self.ops[i][0] == pysam.CINS:
                        mdOffset += self.ops[i][1]
                    i += 1
                self.md.append((refBase, mdOffset))
                mdOffset += 1

    def _getMD(self) -> tuple:
        if self.md == None:
            self._buildMD()
        if self.md == None or self.mdI >= len(self.md):
            return (None, None)
        while self.md[self.mdI][1] < self.opPos:
            self.mdI += 1
            if self.mdI >= len(self.md):
                return (None, None)
        return self.md[self.mdI]

    @property
    def opLength(self) -> int:
        return self.ops[self.opsI][1] if len(self.ops) else 0

    @property
    def opEnd(self) -> int:
        l = self.opLength
        if l == 0:
            return self.opStart
        else:
            return self.opStart + l - 1

    def skipClipped(self, hardOnly: bool = False) -> int:
        count = 0
        if self.opLength > 0:
            while self.opsI < len(self.ops) and (self.op == pysam.CHARD_CLIP if hardOnly else self.clipped):
                count += self.opLength
                self.opsI += 1
            self.opStart = count
            self.opPos = self.opStart

        return count

    def skipToPos(self, pos: int): # Pos is in cigar space
        if self.opPos < 0:
            self.opPos = 0
        dist = pos - self.opPos

        while dist > 0:
            if pos < self.opEnd:
                if self.inSeq:
                    self.seqPos += dist
                if self.inRef:
                    self.refPos += dist
                dist = 0
                self.opPos = pos
            else:
                d = self.opEnd - self.opPos
                dist -= d
                if self.inSeq:
                    self.seqPos += d
                if self.inRef:
                    self.refPos += d
                self.opStart += self.opLength
                self.opPos = self.opStart
                self.opsI += 1
                if self.opsI >= len(self.ops):
                    return False
        return True

    def skipToRefPos(self, pos): # Pos is in reference space
        if self.opPos < 0:
            self.opPos = 0
        dist = pos - self.refPos
        while dist > 0:
            if self.inRef:
                if dist + self.opPos <= self.opEnd +1:
                    if self.inSeq:
                        self.seqPos += dist
                    self.refPos += dist
                    self.opPos += dist
                    dist = 0
                else:
                    dist -= self.opEnd - self.opPos
                    if not self.nextOp(): return False
            else:
                if not self.nextOp(): return False
        return True


    def skipToNonRef(self) -> bool: # Move iterator to next non-reference cigar position (variant in MD tag)
        md = self._getMD()
        if md[0] is None:
            return False
        if md[1] == self.opPos:
            self.mdI += 1
        md = self._getMD()
        if md[0] is None or not self.valid:
            return False
        return self.skipToRefPos(md[1])

    @property
    def inRef(self) -> bool: # Returns true if the passed operation has a reference coordinate
        return self.op in (pysam.CMATCH, pysam.CDEL, pysam.CREF_SKIP, pysam.CEQUAL, pysam.CDIFF)

    @property
    def inSeq(self) -> bool: # Returns true if the passed operation has a sequence coordinate
        return self.op in (pysam.CMATCH, pysam.CINS, pysam.CSOFT_CLIP, pysam.CEQUAL, pysam.CDIFF)

    @property
    def clipped(self):
        return self.op in (pysam.CHARD_CLIP, pysam.CSOFT_CLIP)

    @property
    def refBase(self) -> str:
        return (self.seqBase if self.matchesRef else self._getMD()[0]) if self.inRef else ""

    @property
    def matchesRef(self) -> bool:
        return self._getMD()[1] != self.opPos #self.getSeqBase() == self.getRefBase()

    @property
    def seqBase(self) -> str:
        return self.record.query_sequence[self.seqPos] if self.inSeq else ""

    @seqBase.setter
    def setSeqBase(self, str) -> bool:
        if self.inSeq:
            self.record.query_sequence[self.seqPos] = str
        else:
            return False
        return True

    @property
    def baseQual(self) -> int:
        return self.record.query_qualities[self.seqPos] if self.inSeq else None

    @baseQual.setter
    def setBaseQual(self, qual: int) -> bool:
        if self.inSeq:
            self.record.query_qualities[self.seqPos] = qual
        else:
            return False
        return True

    @property
    def op(self) -> int:
        return self.ops[self.opsI][0]

    @property
    def opRange(self):
        return self.ops[self.opsI]