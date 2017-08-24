import pysam
from sortedcontainers import SortedDict, SortedListWithKey

def BEDtoRef(bedFile: pysam.TabixFile, faiFile: pysam.TabixFile) -> SortedListWithKey:
    index = SortedDict()
    list = SortedListWithKey(key=lambda x: x[0])
    for region in pysam.tabix_iterator(faiFile, pysam.asTuple):
        index[region[0]] = region[2]
    for region in  pysam.tabix_iterator(bedFile, pysam.asBed):
        start = index[region.contig] + region.start
        end = index[region.contig] + region.end
        list.add((start, end))

def inCoords(coord: int, coords: SortedListWithKey) -> bool:
    i = coords.bisect(coord)
    return (coords[i][0] <= coord <= coords[i][1]) \
           or (len(coords) > i and coords[i+1][0] <= coord <= coords[i+1][1]) \
           or (i > 0 and coords[i-1][0] <= coord <= coords[i-1][1])
