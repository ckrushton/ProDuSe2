#!/usr/bin/env python3
import io, time
from sys import stderr, stdout

# If running directly, this works fine
try:
    import FastqRecord
    from configutator import ConfigMap, ArgMap, loadConfig, InvalidParamter
# If installed
except ImportError:
    from ProDuSe import FastqRecord
    from ProDuSe.configutator import ConfigMap, ArgMap, loadConfig  

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

@ConfigMap()
@ArgMap()
def trim(input: list, output: io.IOBase, barcode_distance: int, barcode_sequence: str, reverse: bool = False, verbose: bool = False, logStream: io.IOBase=stderr) -> (int, int):
    """
    Trims barcodes from reads in a fastq file
    :param input: A file or stream handle to read input data
    :param output: A file or stream handle to output data
    :param barcode_distance: The maximum number of differences from the barcode sequence before the read is rejected. Set to negative to output rejected reads only.
    :param barcode_sequence: The IUPAC barcode sequence to match against.
    :param reverse: Set to true to read sequences in reverse, looking for the barcode at the other end.
    :param verbose: Provide verbose output while processing
    :return:
    """

    # Used to print out update messages:
    printPrefix = "PRODUSE-TRIM"
    stdout.write("\t".join([printPrefix, time.strftime('%X'), "Starting...\n"]))

    # Prepare the input files
    mainInput = input[0]
    if len(input) == 1:  # Single-end sequencing or read
        mated = False
    else:  # Paired-end sequencing
        mated = True
        mateInput = input[1]
    if verbose:
        logStream.write("\n") # This will be deleted by the next write

    record1 = FastqRecord.FastqRecord()
    if mated:
        record2 = FastqRecord.FastqRecord()
    invert = False
    count = 0
    discard = 0
    warned = False
    if barcode_distance < 0:
        invert = True
        barcode_distance *= -1
    while record1.read(mainInput) and (not mated or record2.read(mateInput)):
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
        for record, barcode in ((record1, barcode1), (record2, barcode2)) if mated else ((record1, barcode1),):
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
            record.write(output)
        if count % 100000 == 0:
            # Do a quick sanity check: If a significant fraction of reads have been discarded, it is likely that the barcodes is incorrect
            # In this case, inform the user, but continue anyways
            discardRate = float(discard)/float(count)*100
            if not warned and discardRate > 60:
                stdout.write("WARNING: A significant fraction of reads are being discarded. Are you sure the barcode sequence is correct?\n")
                stdout.write("If you are unsure, try running \"produse adapter_predict\" or \"adapter_predict.py\" to double-check the barcode sequence\n")
                warned = True
            if not mated:
                stdout.write("\t".join([printPrefix, time.strftime('%X'), "Reads Processed: " + str(count), "Discard Rate: " + str(discardRate) + "%\n"]))
            else:
                stdout.write("\t".join([printPrefix, time.strftime('%X'), "Pairs Processed: " + str(count), "Discard Rate: " + str(discardRate) + "%\n"]))
    # Now that all the reads have been processed, generate some final stats
    if verbose:
        logStream.write("\x1b[F\x1b[2K\r{file}\tTotal records: {total}\tRecords discarded: {discard}\n".format(file=output.name if hasattr(output, 'name') else 'Streaming', total=count, discard=discard))
    if not mated:
        stdout.write("\t".join([printPrefix, time.strftime('%X'), "Total Reads Processed: " + str(count), "Discard Rate: " + str(float(discard)/float(count)*100) + "%\n"]))
    else:
        stdout.write("\t".join([printPrefix, time.strftime('%X'), "Total Pairs Processed: " + str(count), "Discard Rate: " + str(float(discard)/float(count)*100) + "%\n"]))
    stdout.write("\t".join([printPrefix, time.strftime('%X'), "Trim Complete\n"]))
    output.close()
    return discard, count


def main(args=None):
    from sys import stdin, argv
    import gzip, os, errno
    if args is None:
        args = argv
    """
    pathsDoc = [
        ('reads.fastq[.gz]', 'Fastq with reads to trim'),
        ('[mates.fastq[.gz]]', 'Mates fastq that will also be trimmed and interlaced in the output'),
        ('output.fastq[.gz]', 'Path to output trimmed fastq to')
    ]
    """
    for argmap, paths in loadConfig(args, (trim,), title='Trim V1.0'):

        # Open the input and output files/streams
        inFiles = argmap[trim]["input"]
        # Handle input files
        # Only a maximum of two paired end fastq files can be specified
        if len(inFiles) > 2:
            stderr.write("ERROR: A maximum of two input FASTQ files can be specified\n")
            exit(1)
        argmap[trim]["input"] = []
        for file in inFiles:
            if file == "-":
                # Since a pipe can only be used to provide a single input, error out if both arguments specified pipes
                if stdin in argmap[trim]["input"]:
                    stderr.write("ERROR: A pipe can only specify a single input, not both\n")
                    exit(1)
                argmap[trim]["input"].append(stdin)
            else:
                x = open(file, 'rb')
                if x.peek(2)[:2] == b'\037\213':
                    argmap[trim]["input"].append(gzip.open(file))
                else:
                    argmap[trim]["input"].append(x)

        # Open output file
        outStream = argmap[trim]["output"]
        if outStream == "-":
            argmap[trim]["output"] = stdout
        else:
            if os.path.dirname(outStream) and not os.path.exists(os.path.dirname(outStream)):
                try:
                    os.makedirs(os.path.dirname(outStream))
                except OSError as e:
                    if e.errno != errno.EEXIST:
                        raise
            if outStream.endswith('.gz'):
                argmap[trim]["output"] = gzip.open(outStream, 'wb+')
            else:
                argmap[trim]["output"] = open(outStream, 'wb+')
        trim(**argmap[trim])


if __name__ == "__main__":
    main()
