import pysam
import sortedcontainers
from configutator import ConfigMap

import Levenshtein
from FamilyRecord import FamilyRecord

@ConfigMap(inFile=None, outFile=None)
def collapse(inFile: pysam.AlignmentFile, outFile: pysam.AlignmentFile, barcode_distance: int, barcode_mask: str):
    records = sortedcontainers.SortedDict() #SortedListWithKey(key=lambda x: x.name)
    startPos = 0
    for record in inFile.fetch(until_eof=True):
        if startPos != record.reference_start: # Families are partially determined by start coord
            while len(records):
                family, familyRecord = records.popitem()
                toDel = [] # Can't modify while iterating so store index
                # Search for records that are within threshold hamming distance
                for i in range(len(records)):
                    f, fRec = records.items()[i]
                    if Levenshtein.hamming(family, f) < barcode_distance: # records with family names with a hamming distance < threshold are aggregated
                        toDel.append(i)
                        familyRecord += fRec

                # Remove records once merged
                for i in reversed(toDel): # Delete in reverse to preserve indicies
                    del records.items()[i]

                outFile.write(familyRecord.toPysam())
            startPos = record.reference_start

        # Build family name
        family = "R" if record.is_reverse else "F" # Keep forward and reverse reads separate
        bc = record.get_tag("BC")
        for i in range(len(barcode_mask)):
            if barcode_mask[i] == "1":
                family += bc[i]

        # Add or update family to records
        familyRecord = records.get(family)
        if familyRecord is None:
            familyRecord = FamilyRecord(family, startPos)
            records[family] = familyRecord
        else:
            familyRecord.aggregate(record)

if __name__ == "__main__":
    import getopt
    from sys import stdout, stdin, stderr, argv

    threads = 1
    maxFSize = 0
    outFormat = 'w'
    ordered = False
    verbose = False

    # Command line options
    stderr.write("Collapse v1.0\n")
    ops, paths = getopt.gnu_getopt(argv[1:], 't:m:ho:sv')
    if not len(ops) and not len(paths):
        stderr.write("Use -h for help.\n")
    else:
        for op, val in ops:
            if op == '-h':
                stderr.write("Clip overlapping reads from SAM/BAM/CRAM file\n"
                             "Use: clip.py [-tmos] [input file path | < infile > outfile] [output file path]\n"
                             "If no paths are given stdin and stdout are used.\n"
                             "-t # Threads to use for processing (Default=1)\n"
                             "-m # Maximum family depth (0 for no maximum, Default=0)\n"
                             "-o [sbuc] Output format: s=SAM (Default), b=BAM compressed, bu=BAM uncompressed, c=CRAM\n"
                             "-s Output in coordinate sorted order (Unsorted input may fill RAM), if not set will output in arbitrary order (Minimal RAM)\n"
                             "-v Verbose status output\n")
                exit()
            elif op == '-t':
                threads = int(val or 1)
            elif op == '-m':
                maxFSize = int(val or 0)
            elif op == '-o' and val != 's':
                outFormat += val
            elif op == '-s':
                ordered = True
            elif op == '-v':
                verbose = True