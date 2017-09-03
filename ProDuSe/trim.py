#!/usr/bin/env python3
import io
import gzip
from sys import stdin, stdout, stderr
from configutator import ConfigMap, ArgMap, PositionalArg, TransformCfg
import time
import os
import errno

# If called directly as a script
try:
    import FastqRecord
# If installed
except ImportError:
    from ProDuSe import FastqRecord

IUPACCodeDict = {
    'A': 'A',  # Adenine
    'C': 'C',  # Cytosine
    'G': 'G',  # Guanine
    'T': 'TU',  # Thymine
    'U': 'TU',  # Uracil
    'R': 'AG',  # A or G
    'Y': 'CT',  # C or T
    'S': 'GC',  # G or C
    'W': 'AT',  # A or T
    'K': 'GT',  # G or T
    'M': 'AC',  # A or C
    'B': 'CGT',  # C or G or T
    'D': 'AGT',  # A or G or T
    'H': 'ACT',  # A or C or T
    'V': 'ACG',  # A or C or G
    'N': 'ACGTURYSWKMBDHVN',  # any base
    '.': '.-',  # gap
    '-': '.-'  # gap
}


def IUPACMatch(code1: str, code2: str) -> bool:
    """
    Compares two IUPAC nucleotide codes and returns if code2 is within the IUPAC specified range of code1
    :param code1: IUPAC nucleotide code
    :param code2: IUPAC nucleotide code
    :return: True if the codes are within the specified range, False otherwise
    """
    return code2 in IUPACCodeDict[code1]


def openGZ(path: str, mode: str):
    """
    Open a possibly gzipped file. Checks magic bytes at beginning of file to determine compression.
    :param path: Path to file to open
    :param mode: Mode to open file in (same as open() paramater)
    :return: A open file object, similar to open()
    """
    fh = open(path, mode)
    if fh.peek(2)[:2] == b'\037\213':
        fh = gzip.open(fh)
    return fh


def createAndOpen(path: str, mode: str):
    """
    Open a file, creating the full directory path to the file if necessary.
    :param path: Path to file to open
    :param mode: Mode to open file in (same as open() paramater)
    :return: A open file object, similar to open()
    """
    if os.path.dirname(path) and not os.path.exists(os.path.dirname(path)):
        try:
            os.makedirs(os.path.dirname(path))
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
    if path.endswith('.gz'):
        return gzip.open(path, mode)
    else:
        return open(path, mode)


@ConfigMap(inStream=TransformCfg('input', (openGZ, None, ('rb',))), mateStream=TransformCfg('mate', (openGZ, None, ('rb',))), outStream=TransformCfg('output', (createAndOpen, None, ('wb+',))), logStream=TransformCfg('log', (createAndOpen, None, ('w+',))))
@ArgMap(inStream=PositionalArg(0, 'reads.fastq[.gz]', 'Path to fastq with reads to trim.', (openGZ, None, ('rb+',))), mateStream=PositionalArg(1, 'mates.fastq[.gz]', 'Mates fastq that will also be trimmed and interlaced in the output.', (openGZ, None, ('rb+',))), outStream=PositionalArg(2, 'output.fastq[.gz]', 'Path to output trimmed fastq to.', (createAndOpen, None, ('wb+',))), logStream=PositionalArg(3, 'log', 'Path to output verbose log.', (createAndOpen, None, ('wb+',))))
def trim(barcode_distance: int, barcode_sequence: str, reverse: bool = False, inStream: io.IOBase = stdin, outStream: io.IOBase = stdout, mateStream: io.IOBase = None, verbose: bool = True, logStream: io.IOBase=stderr) -> (int, int):
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

    printPrefix = "PRODUSE-TRIM"
    if verbose:
        logStream.write("\t".join([printPrefix, time.strftime('%X'), "Starting...\n"]))
    record1 = FastqRecord.FastqRecord()
    mated = mateStream is not None
    if mated:
        record2 = FastqRecord.FastqRecord()
    invert = False
    count = 0  # Record counter
    discard = 0  # Discarded records counter
    if barcode_distance < 0:
        invert = True
        barcode_distance *= -1

    while record1.read(inStream) and (not mated or record2.read(mateStream)):
        """
        if verbose:
            logStream.write("\x1b[F\x1b[2K\r{file}\tWorking on record: {record}\tRecords processed: {total}\n".format(file=outStream.name if hasattr(outStream, 'name') else 'Streaming', record=record1.name, total=count))
        """
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
            if invert != (mismatch > barcode_distance):  # XOR
                discard += 1
                # Keep empty record to make BWA happy about having read mate
                record.seq = ''
                record.qual = ''
            else:
                record.trim(0 if reverse else len(barcode_sequence), len(barcode_sequence) if reverse else 0)
            record.write(outStream)

            # Provide status updates
        if verbose and count % 1000 == 0:
            if mated:
                logStream.write("\t".join([printPrefix, time.strftime('%X'), "Read Pairs Trimmed: %s Discard Rate: %s \r" % (count, discard / count)]))
            else:
                logStream.write("\t".join([printPrefix, time.strftime('%X'), "Reads Trimmed: %s Discard Rate: %s \r" % (count, discard / count)]))

    if verbose:
        logStream.write("\t".join(["\n" + printPrefix, time.strftime('%X'), "Trimming Complete\n"]))
    outStream.close()
    return discard, count


if __name__ == "__main__":
    from sys import stdout, stdin, argv
    from configutator import loadConfig

    cfgs = loadConfig(argv, (trim,), title="Trim V1.0")
    while True:  # Skip missing inputs
        try:
            argmap = next(cfgs)
            trim(**argmap[trim])
        except ValueError:
            continue  # TODO add verbose output
        except StopIteration:
            break

