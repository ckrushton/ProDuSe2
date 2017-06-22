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

if __name__ == "__main__":
    pass #TODO