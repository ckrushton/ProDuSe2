#!/usr/bin/env python3
import io, gzip
from sys import stdin, stdout, stderr
import FastqRecord
from configutator import ConfigMap, ArgMap, PositionalArg, TransformCfg

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

def openGZ(path, mode):
    fh = open(path, mode)
    if fh.peek(2)[:2] == b'\037\213':
        fh = gzip.open(fh)
    return fh

def createAndOpen(path, mode):
    if not os.path.exists(os.path.dirname(path)):
        try:
            os.makedirs(os.path.dirname(path))
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
    if path.endswith('.gz'):
        return gzip.open(path, mode)
    else:
        return open(path, mode)

@ConfigMap(inStream=TransformCfg('input', (openGZ, None, ('rb',))), mateStream= TransformCfg('mate', (openGZ, None, ('rb',))), outStream=TransformCfg('output', (createAndOpen, None, ('wb+',))), logStream=TransformCfg('log', (createAndOpen, None, ('w+',))))
@ArgMap(inStream=PositionalArg(0, 'reads.fastq[.gz]', 'Path to fastq with reads to trim.', (openGZ, None, ('wb+',))), mateStream=PositionalArg(1, 'mates.fastq[.gz]', 'Mates fastq that will also be trimmed and interlaced in the output.', (openGZ, None, ('wb+',))), outStream=PositionalArg(2, 'output.fastq[.gz]', 'Path to output trimmed fastq to.', (createAndOpen, None, ('wb+',))), logStream=PositionalArg(3, 'log', 'Path to output verbose log.', (createAndOpen, None, ('wb+',))))
def trim(barcode_distance: int, barcode_sequence: str, reverse: bool = False, inStream: io.IOBase = stdin, outStream: io.IOBase = stdout, mateStream: io.IOBase = None, verbose: bool = False, logStream: io.IOBase=stderr) -> (int, int):
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
            logStream.write("\x1b[F\x1b[2K\r{file}\tWorking on record: {record}\tRecords processed: {total}\n".format(file=outStream.name if hasattr(outStream, 'name') else 'Streaming', record=record1.name, total=count))
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
        logStream.write("\x1b[F\x1b[2K\r{file}\tTotal records: {total}\tRecords discarded: {discard}\n".format(file=outStream.name if hasattr(outStream, 'name') else 'Streaming', total=count, discard=discard))
    outStream.close()
    return discard, count

if __name__ == "__main__":
    from sys import stdout, stdin, argv
    import os, errno
    from configutator import loadConfig
    for argmap in loadConfig(argv, (trim,), title='Trim V1.0'):
        trim(**argmap[trim])
