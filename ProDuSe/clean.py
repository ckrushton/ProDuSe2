import pysam
from io import IOBase

from condense import condense
from clip import mergeRecord, trimRecord

def clean(inStream: IOBase, outStream: IOBase):
    inFile = pysam.AlignmentFile(inStream, 'r', check_sq=False)
    outFile = pysam.AlignmentFile(outStream, "w", template=inFile)
    firstRecord = next(inFile)
    secondRecord = next(inFile)
    try:
        while firstRecord and secondRecord:
            if firstRecord.query_name != secondRecord.query_name:
                raise ValueError("Non-contiguous read mate pair found, exiting.")
            #if regions and (not BEDtoRef.inCoords(firstRecord.reference_start, regions) \
            #        or not BEDtoRef.inCoords(firstRecord.reference_end, regions) \
            #        or not BEDtoRef.inCoords(secondRecord.reference_start, regions) \
            #        or not BEDtoRef.inCoords(secondRecord.reference_end, regions)):
            #    continue

            # Condense BWA output
            condense(firstRecord)
            condense(secondRecord)

            forwardRecord = firstRecord if not firstRecord.is_reverse else secondRecord
            reverseRecord = firstRecord if firstRecord.is_reverse else secondRecord

            mergeRecord(reverseRecord, forwardRecord)
            trimRecord(reverseRecord, forwardRecord.reference_end)

            outFile.write(firstRecord)
            outFile.write(secondRecord)
            firstRecord = next(inFile)
            secondRecord = next(inFile)
    except:
        pass
    finally:
        inFile.close()
        outFile.close()

if __name__ == "__main__":
    import sys
    clean(sys.stdin, sys.stdout)