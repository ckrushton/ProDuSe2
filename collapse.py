import pysam
import sortedcontainers
from configutator import ConfigMap

import Levenshtein
from FamilyRecord import FamilyRecord

@ConfigMap(inFilePath=None, outFilePath=None)
def collapse(barcode_distance: int, barcode_mask: str, inFilePath: str = None, outFilePath: str = None):
    inFile = pysam.AlignmentFile(inFilePath or sys.stdin)
    outFile = pysam.AlignmentFile(outFilePath or sys.stdout, "rb")
    records = sortedcontainers.SortedDict() #SortedListWithKey(key=lambda x: x.name)
    startPos = 0
    for record in inFilePath.fetch(until_eof=True):
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

                outFilePath.write(familyRecord.toPysam())
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
    from configutator import handleArgs
    import sys
    for argmap in handleArgs(sys.argv, (collapse,)):
        global config
        config = argmap
        collapse(**argmap[collapse])
        break