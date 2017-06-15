import pysam
import sortedcontainers

import Levenshtein
from FamilyRecord import FamilyRecord

def collapse(inFile: pysam.AlignmentFile, outFile: pysam.AlignmentFile, threshold: int, mask: str):
    records = sortedcontainers.SortedListWithKey(key=lambda x: x.name)
    startPos = 0
    for record in inFile.fetch():
        if startPos != record.reference_start: # Families are partially determined by start coord
            while len(records):
                family, familyRecord = records.popitem()
                toDel = [] # Can't modify while iterating so store index
                # Search for records that are within threshold hamming distance
                for i in range(len(records)):
                    f, fRec = records.items()[i]
                    if Levenshtein.hamming(family, f) < threshold: # records with family names with a hamming distance < threshold are aggregated
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
        for i in range(len(mask)):
            if mask[i] == "1":
                family += bc[i]

        # Add or update family to records
        familyRecord = records.get(family)
        if familyRecord is None:
            familyRecord = FamilyRecord(family, startPos)
            records.add(familyRecord)
        else:
            familyRecord.aggregate(record)

def getArgs(parser):
    import configargparse
    # Collapse arguments
    collapseArgs = parser.add_argument_group("Collapse Arguments")

    collapseArgs.add(
        "-dp", "--duplex_position",
        metavar="STR",
        type=str,
        required=True,
        help="The positions in the adapter sequence to include in distance calculations between forward and reverse reads, 0 for no, 1 for yes"
    )

    # Used to maintain backwards compatibility with the poorly-named strand mis-match
    adapterMismatch = collapseArgs.add_mutually_exclusive_group(required=True)
    adapterMismatch.add(
        "-amm", "--adapter_max_mismatch",
        type=int,
        help="The maximum number of mismatches allowed between the expected and actual adapter sequences",
    )
    adapterMismatch.add(
        "--strand_max_mismatch",
        type=int,
        help=configargparse.SUPPRESS,
    )

    collapseArgs.add(
        "-dmm", "--duplex_max_mismatch",
        type=int,
        required=True,
        help="The maximum number of mismatches allowed between the expected and actual duplexed adapter sequences",
    )

    collapseArgs.add(
        "-smm", "--sequence_max_mismatch",
        type=int,
        required=False,
        default=20,
        help="The maximum number of mismatches allowed in an alignment"
    )

    collapseArgs.add(
        "-oo", "--original_output",
        type=str,
        required=False,
        action="append",
        default=None,
        help="A pair of empty fastq files to rewrite original fastqs with duplex information"
    )

    collapseArgs.add(
        "-sf", "--stats_file",
        type=str,
        required=False,
        default=None,
        help="An optional output file to list stats generated during collapse"
    )

