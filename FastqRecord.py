import io

class FastqRecord(object):
    "Parses an input stream containing well formed fastq records"

    # Name Description1 Sequence Description2 Quality
    name, desc1, seq, desc2, qual = "","","","",""
    
    def read(self, stream:io.IOBase) -> bool:
        line1 = stream.readline()
        if line1 is None: return False
        line2 = stream.readline()
        if line2 is None: return False
        line3 = stream.readline()
        if line3 is None: return False
        line4 = stream.readline()
        if line4 is None: return False

        # Generate index for record data
        nameEnd = self.line1.find(' ')-1
        self.name = line1[1:nameEnd]
        self.desc1 = line1[nameEnd+2:-1]
        self.seq = line2[:-1]
        self.desc2 = line3[self.line3.find(' ')+1:-1]
        self.qual = line4[:-1]
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
        self.seq = self.seq[left:-right]
        self.qual = self.qual[left:-right]

    def write(self, stream: io.IOBase):
        stream.writelines(['@', self.name, ' ' + self.desc1 if self.desc1 != '' else '', '\n',
                          self.seq, '\n',
                           '+ ' + self.desc2 if self.desc2 != '' else '+', '\n',
                           self.qual, '\n'])