#!/usr/bin/env python3
import io
from sys import stderr
import FastqRecord
from configutator import ConfigMap, ArgMap

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

@ConfigMap(inStream=None, mateStream=None, outStream=None, logStream=None)
@ArgMap(inStream=None, mateStream=None, outStream=None, logStream=None)
def trim(inStream: io.IOBase, outStream: io.IOBase, barcode_distance: int, barcode_sequence: str, reverse: bool = False, mateStream: io.IOBase = None, verbose: bool = False, logStream: io.IOBase=stderr) -> (int, int):
    """
    Trims barcodes from reads in a fastq file
    :param inStream: A file or stream handle to read input data
    :param mateStream: A file or stream handle to read mate input data
    :param outStream: A file or stream handle to output data
    :param barcode_distance: The maximum number of differences from the barcode sequence before the read is rejected. Set to negative to output rejected reads only.
    :param barcode_sequence: The IUPAC barcode sequence to match against.
    :param reverse: Set to true to read sequences in reverse, looking for the barcode at the other end.
    :param verbose: Provide verbose output while processing
    :return:
    """
    if verbose:
        logStream.write("\n") # This will be deleted by the next write
    record1 = FastqRecord.FastqRecord()
    mated = mateStream is not None
    if mated:
        record2 = FastqRecord.FastqRecord()
    invert = False
    count = 0
    discard = 0
    if barcode_distance < 0:
        invert = True
        barcode_distance *= -1
    while record1.read(inStream) and (not mated or record2.read(mateStream)):
        if verbose:
            logStream.write("\x1b[F\x1b[2K\r{file}\tWorking on record: {record}\tRecords processed: {total}\n".format(file=inStream.name if hasattr(inStream, 'name') else 'Streaming', record=record1.name, total=count))
        count += 1
        # Store barcode in sequence description
        barcode1 = (record1.seq[-len(barcode_sequence):] if reverse else record1.seq[:len(barcode_sequence)])
        if mated:
            barcode2 = (record2.seq[-len(barcode_sequence):] if reverse else record2.seq[:len(barcode_sequence)])
        else:
            barcode2 = ''
        record1.desc1 += '\tBC:Z:' + barcode1 + barcode2
        if mated:
            record2.desc1 += '\tBC:Z:' + barcode2 + barcode1
        for record, barcode in ((record1, barcode1), (record2, barcode2)) if mated else ((record1, barcode1)):
            mismatch = 0
            mismatch += sum(not IUPACMatch(bc, sq) for bc, sq in zip(barcode_sequence, reversed(barcode) if reverse else barcode))
            if len(barcode) < len(barcode_sequence):
                mismatch += len(barcode_sequence) - len(barcode)
            if invert != (mismatch > barcode_distance): #XOR
                discard += 1
                # Keep empty record to make BWA happy about having read mate
                record.seq = ''
                record.qual = ''
            else:
                record.trim(0 if reverse else len(barcode_sequence), len(barcode_sequence) if reverse else 0)
            record.write(outStream)
    if verbose:
        logStream.write("\x1b[F\x1b[2K\r{file}\tTotal records: {total}\tRecords discarded: {discard}\n".format(file=inStream.name if hasattr(inStream, 'name') else 'Streaming', total=count, discard=discard))
    outStream.close()
    return discard, count

if __name__ == "__main__":
    from sys import stdout, stdin, argv
    import gzip, os, errno
    from configutator import loadConfig
    pathsDoc = [
        ('reads.fastq[.gz]', 'Fastq with reads to trim'),
        ('[mates.fastq[.gz]]', 'Mates fastq that will also be trimmed and interlaced in the output'),
        ('output.fastq[.gz]', 'Path to output trimmed fastq to')
    ]
    for argmap, paths in loadConfig(argv, (trim,), title='Trim V1.0', positionalDoc=pathsDoc):
        if len(paths):
            inFile = open(paths[0], 'rb')
            if inFile.peek(2)[:2] == b'\037\213':
                inFile = gzip.open(inFile)
        else:
            inFile = stdin
        if len(paths) > 2:
            mateFile = open(paths[1], 'rb')
            if mateFile.peek(2)[:2] == b'\037\213':
                mateFile = gzip.open(mateFile)
        else:
            mateFile = None
        if len(paths) > 1:
            i = 2 if len(paths) > 2 else 1
            if not os.path.exists(os.path.dirname(paths[i])):
                try:
                    os.makedirs(os.path.dirname(paths[i]))
                except OSError as e:
                    if e.errno != errno.EEXIST:
                        raise
            if paths[i].endswith('.gz'):
                outFile = gzip.open(paths[i], 'wb+')
            else:
                outFile = open(paths[i], 'wb+')
        else:
            outFile = stdout
        trim(inFile, outFile, mateStream=mateFile, **argmap[trim])
