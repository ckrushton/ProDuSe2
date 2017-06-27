import io
import FastqRecord
from configutator import ConfigMap

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

@ConfigMap(inStream=None, outStream=None)
def trim(inStream: io.IOBase, outStream: io.IOBase, barcode_distance: int, barcode_sequence: str, reverse: bool = False) -> (int, int):
    record = FastqRecord.FastqRecord()
    invert = False
    count = 0
    discard = 0
    if barcode_distance < 0:
        invert = True
        barcode_distance *= -1
    while record.read(inStream):
        mismatch = 0
        count += 1
        # Store barcode in sequence description
        record.desc1 += '\tBC:Z:' + (record.seq[-len(barcode_sequence):] if reverse else record.seq[:len(barcode_sequence)])
        for i in range(0, len(barcode_sequence), 1):
            mismatch += not IUPACMatch(barcode_sequence[i], record.seq[-i if reverse else i])
        if invert != mismatch > barcode_distance: #XOR
            discard += 1
            # Keep empty record to make BWA happy about having read mate
            record.seq = ''
            record.qual = ''
        record.write(outStream)
    return discard, count

if __name__ == "__main__":
    pass #TODO