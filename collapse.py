import pysam
import sortedcontainers
from merge import mergeRecord
import Levenshtein

def collapse(inFile: pysam.AlignmentFile, outFile: pysam.AlignmentFile, threshold: int, mask: str):
    records = sortedcontainers.SortedDict()
    startPos = 0
    for record in inFile.fetch():
        if startPos != record.reference_start: # Families are partially determined by start coord
            while len(records):
                family, familyRecord = records.popitem()
                toDel = [] # Can't modify while iterating so store index
                for i in range(len(records)):
                    f, r = records.items()[i]
                    if Levenshtein.hamming(family, r[0]) < threshold: # records with family names with a hamming distance < threshold are merged
                        toDel.append(i)
                        mergeRecord(r[0], familyRecord[0])
                        familyRecord[1] += r[1]
                        # TODO copy mapping Q?
                for i in reversed(toDel): # Delete in reverse to preserve indicies
                    del records.items()[i]
                familyRecord[0].query_name = family
                outFile.write(familyRecord[0])
            startPos = record.reference_start
        family = "R" if record.is_reverse else "F" # Keep forward and reverse reads separate
        bc = record.get_tag("BC")
        for i in range(len(mask)):
            if mask[i] == "1":
                family += bc[i]

        familyRecord = records.get(family)
        if familyRecord is None:
            records[family] = [record, 1]
        else:
            familyRecord[1] += 1
            mergeRecord(record, familyRecord[0])
