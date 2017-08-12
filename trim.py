#!/usr/bin/env python3
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
def trim(inStream: io.IOBase, outStream: io.IOBase, barcode_distance: int, barcode_sequence: str, reverse: bool = False, verbose: bool = False) -> (int, int):
    """
    Trims barcodes
    :param inStream: A file or stream handle to read input data
    :param outStream: A file or stream handle to output data
    :param barcode_distance: The maximum number of differences from the barcode sequence before the read is rejected. Set to negative to output rejected reads only.
    :param barcode_sequence: The IUPAC barcode sequence to match against.
    :param reverse: Set to true to read sequences in reverse, looking for the barcode at the other end.
    :param verbose: Provide verbose output while processing
    :return:
    """
    record = FastqRecord.FastqRecord()
    invert = False
    count = 0
    discard = 0
    if barcode_distance < 0:
        invert = True
        barcode_distance *= -1
    while record.read(inStream):
        if verbose:
            stderr.write("\x1b[2K\r{file}\tWorking on record: {record}\tRecords processed: {total}".format(file=outStream.name if hasattr(outStream, 'name') else 'Streaming', record=record.name, total=count))
        mismatch = 0
        count += 1
        # Store barcode in sequence description
        record.desc1 += '\tBC:Z:' + (record.seq[-len(barcode_sequence):] if reverse else record.seq[:len(barcode_sequence)])
        mismatch += sum(not IUPACMatch(bc, sq) for bc, sq in zip(barcode_sequence, reversed(record.seq) if reverse else record.seq))
        if len(record.seq) < len(barcode_sequence):
            mismatch += len(barcode_sequence) - len(record.seq)
        if invert != (mismatch > barcode_distance): #XOR
            discard += 1
            # Keep empty record to make BWA happy about having read mate
            record.seq = ''
            record.qual = ''
        record.write(outStream)
    if verbose:
        stderr.write("\x1b[2K\r{file}\tTotal records: {total}\tRecords discarded: {discard}\n".format(file=outStream.name if hasattr(outStream, 'name') else 'Streaming', total=count, discard=discard))
    return discard, count

if __name__ == "__main__":
    from sys import stdout, stdin, stderr, argv
    import gzip, os, errno
    from configutator import loadConfig
    for argmap, paths in loadConfig(argv, (trim,), title="Trim V1.0"):
        if len(paths):
            inFile = open(paths[0], 'rb')
            if inFile.peek(2)[:2] == b'\037\213':
                inFile = gzip.open(inFile)
        else:
            inFile = stdin
        if len(paths) > 1:
            if not os.path.exists(os.path.dirname(paths[1])):
                try:
                    os.makedirs(os.path.dirname(paths[1]))
                except OSError as e:
                    if e.errno != errno.EEXIST:
                        raise
            if paths[1].endswith('.gz'):
                outFile = gzip.open(paths[1], 'wb+')
            else:
                outFile = open(paths[1], 'wb+')
        else:
            outFile = stdout
        trim(inFile, outFile, **argmap[trim])
