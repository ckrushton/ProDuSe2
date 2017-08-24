import ctypes
class CommaSeparatedList(list):
    def __init__(self, values, value_type=int):
        if isinstance(values, list):
            super().__init__(values)
        elif isinstance(values, str):
            super().__init__([value_type(val) for val in values.split(',')])

    def __str__(self):
        return ",".join(self)

class Tag(ctypes.Structure):
    _fields_ = "TAG", "VALUE"
    def __init__(self, tag):
        if isinstance(tag, str):
            self.TAG, t, self.VALUE = tag.split(':', 3)
            if t in 'AZ':
                pass
            elif t == 'i':
                self.VALUE = int(self.VALUE)
            elif t == 'f':
                self.VALUE = float(self.VALUE)
            elif t == 'H':
                raise NotImplementedError()
            elif t == 'B':
                self.VALUE = CommaSeparatedList(self.VALUE)
        elif isinstance(tag, tuple):
            if len(tag) == 2:
                self.TAG, self.VALUE = tag
                if isinstance(self.VALUE, int):
                    self.TYPE = 'i'
                elif isinstance(self.VALUE, float):
                    self.TYPE = 'f'
                elif isinstance(self.VALUE, list) and len(self.VALUE) and isinstance(self.VALUE[0], int):
                    self.TYPE = 'B'
                elif isinstance(self.VALUE, int):
                    self.TYPE = 'i'
                else:
                    self.TYPE = 'Z'

            else:
                self.TAG, self.VALUE, self.TYPE = tag

class SAMRecord(ctypes.Structure):
    _fields_ = "query_name", "flag", "reference_name", "reference_start", "mapping_quality", "cigartuples", "next_reference_name", "next_reference_start", "template_length", "query_sequence", "query_qualities", "_tags", "_tags_unparsed"
    _slot_parsers = (str, int, str, int, int, list, str, int, int, str, CommaSeparatedList)
    def __init__(self, other: 'SAMRecord' = None):
        if isinstance(other, str):
            currentPos = 0
            for i in range(11):
                nextPos = other.find('\t', currentPos)
                setattr(self, self.__slots__[i], self._slot_parsers[i](other[currentPos + 1:nextPos - 1]))
                currentPos = nextPos
            self._tags_unparsed = other[currentPos+1:]
        elif other:
            for attr in self.__slots__:
                if attr == '_tags' and other.__name__ == "AlignedSegment":
                    setattr(self, attr, other.get_tags(True)) #Clone the passed pysam.AlignedSegment
                else:
                    setattr(self, attr, getattr(other, attr))

    def __str__(self):
        cigar = ""
        qual = ""
        if self.tags is not None:
            self._tags_unparsed = "\t" + "\t".join(["{tag}:{type}:{value}".format(tag, t, val) for tag, val, t in self.tags])
        return "{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}\t{tags}".format(
            self.query_name, self.flag, self.reference_name, self.reference_start+1, self.mapping_quality, cigar, self.next_reference_name,
            self.next_reference_start+1, self.template_length, self.query_sequence, qual) + self.tags_unparsed

    def prepareRewritten(self, other: 'SAMRecord'):
        self.query_name = other.query_name
        self.reference_name = other.reference_name
        self.next_reference_name = other.next_reference_name
        self.next_reference_start = other.next_reference_start