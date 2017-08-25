#!/usr/bin/env python3
import platform
if platform.python_implementation() == 'PyPy':
    import pypysam
else:
    import pysam

import sortedcontainers
import networkx as nx
from configutator import ConfigMap, ArgMap
import io, itertools, collections
from sys import maxsize, stderr

from FamilyRecord import FamilyRecord, appendOrInc

strandThreshold = 0.4 # 40% percent tolerance of deviation for the strand identifier sequence
strandMaskId = 'S'
OSTagName = 'OS' # Id for original start position tag

def getBarcode(record: pysam.AlignedSegment):
    try:
        return record.get_tag("BC")
    except KeyError:
        raise ValueError("Record found without BC tag. All records must have BC tag.")

def pairwise(iterable):
    """Taken from https://docs.python.org/3.4/library/itertools.html#itertools-recipes"""
    #s -> (s0,s1), (s1,s2), (s2, s3), ...
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)

def hamming(str1, str2, max: int = maxsize) -> int:
    """
    Calculate the hamming distance of two strings
    :param str1: The first string to compare
    :param str2: The second string to compare
    :param max: The maximum distance before stopping
    :return: The number of differences between the two strings
    """
    count = 0
    for c1, c2 in zip(str1, str2):
        if c1 != c2: count += 1
        if count >= max: return count
    return count

def maskString(string: str, mask: str, maskCode: str, exclude: bool=False) -> str:
    """
    Return a substring masked by maskCodes
    :param string: The string to mask
    :param mask: A string containing codes to align to string
    :param maskCode: The elements in mask to consider
    :param exclude: If True, returns all elements in string that are not masked by maskCodes, otherwise returns all elements in string where mask macthes maskCode
    :return: A substring of string
    """
    result = ""
    for element, code in zip(string, mask):
        if code in maskCode != exclude: # XOR
            result += element
    return result

def familyName(record: pysam.AlignedSegment, barcodeMask: str) -> str:
    """
    Builds family name from forward/reverse and barcode
    :param record: pysam.AlignedSegment with the BC tag set
    :param barcodeMask: String mask containing '0' for any elements to exclude from the name
    :return: String containing family name
    """
    return ("R" if record.is_reverse else "F") + maskString(getBarcode(record), barcodeMask, '0', True)

def strandCode(familyName: str, mask: str, invert: bool = False) -> str:
    """
    Gets strand specific region of family name
    :param familyName: String result of familyName()
    :param mask: String specifying the location of the strand specific region with 'S'
    :param invert: If True, return family name without strandCode
    :return: String containing the strand specific region of the family name
    """
    mask = mask.replace('0', '') # Family name is result of masking '0' so remove
    return maskString(familyName, mask, strandMaskId, invert)

class Families:
    def __init__(self, barcode_mask: str, barcode_distance: int):
        self._familyRecords = sortedcontainers.SortedDict()
        self._familyMateMap = {}
        self._barcode_mask = barcode_mask
        self._barcode_distance = barcode_distance
        self._recordCounter = 1

    def addRecord(self, record: pysam.AlignedSegment):
        try:
            startPos = record.get_tag(OSTagName)  # type: int
        except KeyError:
            startPos = record.reference_start
        pos = self._familyRecords.get(startPos)
        family = familyName(record, self._barcode_mask)
        if pos is None:
            familyRecord = FamilyRecord(family, startPos, getBarcode(record), record)
            self._familyRecords[startPos] = {family: familyRecord}
        else:
            familyRecord = pos.get(family)
            if familyRecord:
                familyRecord.aggregate(record)
            else:
                familyRecord = FamilyRecord(family, startPos, getBarcode(record), record)
                pos[family] = familyRecord
        self._familyMateMap[(record.query_name, record.is_reverse)] = familyRecord

    def _cutTreeToThreshold(self, forest, tree, thresholdDiameter):
        for branch in nx.connected_component_subgraphs(tree): # tree may already be disconnected
            eccentricity = nx.eccentricity(branch)#, sp=forest.graph['paths'])
            if thresholdDiameter <= nx.diameter(branch, eccentricity):
                peripheral = nx.periphery(branch, eccentricity)
                candidateEdges = {}
                # Determine all candidate edges to remove that will result in a peripheral node being within threshold
                for node1 in peripheral:
                    for node2 in peripheral:
                        if node1 == node2: continue
                        _, path = nx.bidirectional_dijkstra(branch, node1, node2)
                        length = 0
                        for n1, n2 in pairwise(path):
                            data = branch.get_edge_data(n1, n2)
                            length += data['weight']
                            if length > thresholdDiameter:
                                candidateEdges[(n1, n2)] = data
                                break
                maxRankEdge = None
                # Determine candidate edge that will result in maximum total node weight
                for edge, data in candidateEdges.items():
                    branch.remove_edge(*edge)
                    for component in nx.connected_components(branch):
                        rank = 0
                        for node in component:
                            rank += forest.node[node]['weight']
                        if 'rank' not in data or rank > data['rank']:
                            data['rank'] = rank
                    if maxRankEdge is None or maxRankEdge[1] < data['rank']:
                        maxRankEdge = (edge, data['rank'])
                    branch.add_edge(*edge, data)
                if maxRankEdge:
                    forest.remove_edge(*maxRankEdge[0]) #Remove edge from forest to store global result
                    branch.remove_edge(*maxRankEdge[0])
                    # Recurse for each new subtree to ensure it is within threshold
                    for component in nx.connected_component_subgraphs(branch, False):
                        self._cutTreeToThreshold(forest, component, thresholdDiameter)

    def _collapsePosition(self, pos: dict):
        #Collapse families at pos to weighted minimum spanning forest with trees that have no path greater than 2*threshold, maximising family node weight.
        if pos:
            keys = list(pos.keys())
            graph = nx.empty_graph(len(pos)) #type: nx.Graph
            for node, data in graph.nodes_iter(True):
                data['weight'] = pos[keys[node]].size
            # TODO ensure families are seperated by strand threshold and F/R
            for u, v in itertools.combinations(range(len(pos)), 2):
                # Keep families with F/R separate
                if keys[u][0] != keys[v][0]: continue
                # Keep families with different strand specific tags separate
                if hamming(strandCode(keys[u], self._barcode_mask), strandCode(keys[v], self._barcode_mask)) / self._barcode_mask.count(strandMaskId) > strandThreshold: continue
                # Keep families with a hamming distance greater than threshold separate
                dist = hamming(strandCode(keys[u], self._barcode_mask, True), strandCode(keys[v], self._barcode_mask, True), self._barcode_distance)
                if dist < self._barcode_distance: graph.add_edge(u, v, weight=dist)
            thresholdDiameter = 2 * self._barcode_distance
            allShortestPaths = nx.all_pairs_dijkstra_path_length(graph)
            forest = nx.minimum_spanning_tree(graph) #type: nx.Graph # TODO correct greedy algorithm for minimum eccentricity
            forest.graph['paths'] = allShortestPaths
            self._cutTreeToThreshold(forest, forest, thresholdDiameter)
            for component in nx.connected_components(forest):
                component = list(component)
                first = pos[keys[component[0]]]  #type: FamilyRecord
                #Merge together members of component
                for node in component[1:]:
                    other = pos[keys[node]] #type: FamilyRecord
                    for member in other.members: self._familyMateMap[member] = first # Update family mappings
                    first += other
                    del pos[keys[node]]

            #i = 0
            #while len(pos) > i:
            #    for key in keys[i+1:]:
            #        if hamming(key, keys[i], self._barcode_distance) < self._barcode_distance: # records with family names with a hamming distance < threshold are aggregated
            #            pos[keys[i]] += pos[key]
            #            del pos[key]


    def getFamiliesAt(self, pos) -> dict:
        p = self._familyRecords.get(pos)
        if p:
            self._collapsePosition(p)
            return p
        else:
            return {}

    def setMate(self, family: FamilyRecord):
        if family.mate:
            return
        counts = {}
        popular = None
        for member in family.members:
            mateFamily = self._familyMateMap.get(member)
            if mateFamily and not mateFamily.mate:
                if mateFamily not in counts: counts[mateFamily] = 1
                else: counts[mateFamily] += 1
                if not mateFamily.mate and (popular is None or counts[popular] < counts[mateFamily]):
                    popular = mateFamily
        if popular:
            family.mate = popular
            family.mate.mate = family

    def __iter__(self):
        return self._familyRecords.__iter__()

    def __len__(self):
        return len(self._familyRecords)

    def __delitem__(self, pos):
        del self._familyRecords[pos]

    def toPysam(self, family: FamilyRecord) -> (pysam.AlignedSegment, pysam.AlignedSegment):
        self.setMate(family)
        if family.mate.name != family.name:
            family.name = self._recordCounter
            family.mate.name = family.name
            self._recordCounter += 1

        seq = ""
        qual = []
        ops = []
        fQ = []
        fC = []

        # Compile sequence and added tags
        # fQ:B:C Integer array containing Phred score of wrong base chosen during collapse
        # fC:B:I Integer array containing pairs (simmilar to a CIGAR) of count and depth representing, from the start of the alignment, the depth of the family at a position
        for op, base, q, c, q in family:
            appendOrInc(ops, [op, 1])
            if seq:
                seq += base
                qual.append(q)
                fC.append(c)
                fQ.append(q)
        tags = [('BC', family.barcode)]
        if fQ:
            tags.append(('fQ', fQ))
        if fC:
            tags.append(('fC', fC))

        # Copy into pysam record
        record = pysam.AlignedSegment()
        record.query_name = family.name
        record.reference_start = family.pos
        record.mapping_quality = family.maxMapQ
        record.cigartuples = ops
        record.reference_id = family.refID
        record.is_forward = family.forwardCounter >= 0
        if family.mate:
            record.is_paired = True
            record.is_proper_pair = True
            record.next_reference_start = family.mate.pos
            record.next_reference_id = family.mate.refID
            record.mate_is_reverse = family.mate.forwardCounter >= 0
        else:
            record.is_unmapped = True

        record.query_sequence = seq
        record.query_qualities = qual
        if tags:
            record.set_tags(tags)
        return record

def CoordinateSortedInputFamilyIterator(inFile, families):
    lastStartPos = 0
    for record in inFile.fetch(until_eof=True):
        try:
            startPos = record.get_tag(OSTagName)  # type: int
        except KeyError:
            startPos = record.reference_start
        if startPos != lastStartPos:
            f = families.getFamiliesAt(startPos).items()
            for _, family in f:
                yield family
            try:
                if f:
                    del families[startPos]
            except KeyError:
                pass
            lastStartPos = startPos
        families.addRecord(record)

@ConfigMap(inStream=None, outStream=None, logStream=None)
@ArgMap(inStream=None, outStream=None, logStream=None)
def collapse(inStream: io.IOBase, outStream: io.IOBase, barcode_distance: int, barcode_mask: str, verbose: bool=False, logStream: io.IOBase=stderr):
    """
    Collapse reads sharing the same start position, forward/reverse, and barcode into a single family record.
    :param inStream: A file or stream handle to read input data
    :param outStream: A file or stream handle to output data
    :param barcode_distance: The maximum number of differences allowed to be considered the same family
    :param barcode_mask: String mask to denote the positions within a barcode to include in the family name
    :param verbose: Provide verbose output while processing
    :param logStream: The stream to write log output to, defaults to stderr
    :return: None
    """
    if verbose:
        logStream.write("\n") # This will be deleted by the next write
    inFile = pysam.AlignmentFile(inStream)
    outFile = pysam.AlignmentFile(outStream, "wbu" if hasattr(outStream, 'name') else "wb", template=inFile) # compress the data if it is a file
    families = Families(barcode_mask*2, barcode_distance) # reads have mate barcode concatenated so double up mask
    count = 0
    fcount = 0
    lastFCount = 0
    lastFamilyPos = 0
    for family in CoordinateSortedInputFamilyIterator(inFile, families): #type: FamilyRecord
        outFile.write(family.toPysam())
        count += len(family)
        fcount += 1
        if verbose and lastFamilyPos != family.pos:
            logStream.write("\x1b[F\x1b[2K\r{file}\tRecords processed: {count}\t"
                            "Positions buffered: {pos}\tFamily records: {fcount}\t"
                            "Families at last position: {lastFCount}\n".format(
                                file=inStream.name if hasattr(inStream, 'name') else 'Streaming',
                                count=count, pos=len(families),
                                fcount=fcount,
                                lastFCount=fcount-lastFCount)
                            )
            lastFamilyPos = family.pos
            lastFCount = fcount

    if verbose:
        logStream.write("Completed.\n")
    inFile.close()
    outFile.close()
    outStream.close()

if __name__ == "__main__":
    from sys import stdout, stdin, stderr, argv, maxsize, maxsize
    import os, errno
    from configutator import loadConfig
    for argmap, paths in loadConfig(argv, (collapse,), title="Collapse V1.0"):
        if len(paths) and not os.path.exists(paths[0]):
            if 'verbose' in argmap[collapse] and argmap[collapse]['verbose']:
                stderr.write("{} not found, skipping.\n".format(paths[0]))
            continue
        if len(paths) > 1 and not os.path.exists(os.path.dirname(paths[1])):
            try:
                os.makedirs(os.path.dirname(paths[1]))
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise
        collapse(open(paths[0], 'rb') if len(paths) else stdin, open(paths[1], 'wb+') if len(paths) > 1 else stdout, **argmap[collapse])