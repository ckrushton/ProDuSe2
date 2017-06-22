import io
import FastqRecord

IUPACCodeDict = {
    'A' : ['A'],    #Adenine
    'C'	: ['C'],    #Cytosine
    'G'	: ['G'],    #Guanine
    'T' : ['T', 'U'],    #Thymine
    'U' : ['T', 'U'],    #Uracil
    'R' : ['A', 'G'],    #A or G
    'Y' : ['C', 'T'],	#C or T
    'S' : ['G', 'C'],	#G or C
    'W' : ['A', 'T'],	#A or T
    'K' : ['G', 'T'],	#G or T
    'M' : ['A', 'C'],	#A or C
    'B' : ['C', 'G', 'T'],	#C or G or T
    'D' : ['A', 'G', 'T'],	#A or G or T
    'H' : ['A', 'C', 'T'],	#A or C or T
    'V' : ['A', 'C', 'G'],	#A or C or G
    'N' : ['A','C','G','T','U','R','Y','S','W','K','M','B','D','H','V','N'],	#any base
    '.' : ['.', '-'],    #gap
    '-' : ['.', '-']     #gap
}

def IUPACMatch(code1: str, code2: str) -> bool:
    return code2 in IUPACCodeDict[code1]

def trim(inStream: io.IOBase, outStream: io.IOBase, threshold: int, sequence: str, reverse: bool = False) -> (int, int):
    record = FastqRecord.FastqRecord()
    invert = False
    count = 0
    discard = 0
    if threshold < 0:
        invert = True
        threshold *= -1
    while record.read(inStream):
        mismatch = 0
        count += 1
        # Store barcode in sequence description
        record.desc1 += ' BC:Z:' + record.seq[-len(sequence):] if reverse else record.seq[:len(sequence)]
        for i in range(0, len(sequence), 1):
            mismatch += IUPACMatch(sequence[i], record.seq[-i if reverse else i])
        if invert != mismatch > threshold: #XOR
            discard += 1
            # Keep empty record to make BWA happy about having read mate
            record.seq = ''
            record.qual = ''
        record.write(outStream)
    return discard, count

if __name__ == "__main__":
    pass #TODO