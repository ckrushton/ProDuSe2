import io

class FastqRecord(object):
    """Parses an input stream containing well formed fastq records"""
    __slots__ = 'name', 'desc1', 'seq', 'desc2', 'qual'
    def __init__(self):
        # Name Description1 Sequence Description2 Quality
        self.name, self.desc1, self.seq, self.desc2, self.qual = "","","","",""
    
    def read(self, stream:io.IOBase, SAMify:bool = True) -> bool:
        line1 = ' '
        while line1.isspace():  # Ignore white space
            line1 = stream.readline().decode('ascii')
        if line1 == '': return False
        line2 = stream.readline().decode('ascii')
        if line2 == '': return False
        line3 = stream.readline().decode('ascii')
        if line3 == '': return False
        line4 = stream.readline().decode('ascii')
        if line4 == '': return False

        # Generate index for record data
        nameEnd = line1.find(' ')
        self.name = line1[1:nameEnd]
        self.desc1 = line1[nameEnd+1:-1]
        self.seq = line2[:-1]
        nameEnd = line3.find(' ')
        self.desc2 = line3[nameEnd+1:-1] if nameEnd > 0 else ""
        self.qual = line4[:-1]

        if SAMify:
            if self.desc1:
                self.desc1 = "CO:Z:" + self.desc1
            if self.desc2:
                self.desc1 = "CO:Z:" + self.desc2
        return True
    
    def clip(self, start: int, stop: int):
        "Clip a range from the sequence and quality strings"
        if start < 0: start = 0
        if stop >= len(self.seq): stop = len(self.seq)-1
        self.seq = self.seq[:start] + self.seq[stop:]
        self.qual = self.qual[:start] + self.qual[stop:]

    def trim(self, left: int, right: int):
        "Cut left# bases from the beginning of the sequence and quality string, and right# bases from the end"
        if left < 0: left=0
        if right < 0: right=0
        if left >= len(self.seq) or right >= len(self.seq):
            self.seq = ""
            self.qual = ""
            return
        self.seq = self.seq[left:-right or None]
        self.qual = self.qual[left:-right or None]

    def write(self, stream: io.IOBase):
        stream.writelines([b'@', self.name.encode('ascii'), b' ' + self.desc1.encode('ascii') if self.desc1 != '' else '', b'\n',
                          self.seq.encode('ascii'), b'\n',
                           b'+ ' + self.desc2.encode('ascii') if self.desc2 != '' else b'+', b'\n',
                           self.qual.encode('ascii'), b'\n'])