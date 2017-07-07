import pysam
import sortedcontainers

f = pysam.AlignmentFile('/home/ncm3/data/normal.bam')
c = sortedcontainers.SortedListWithKey(key=lambda x: x[1].reference_start)
i = 0
for record in f:
    c.add((i, record))
    i += 1
    print(len(c))
