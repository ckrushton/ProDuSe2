#!/usr/bin/env python3
import platform
if platform.python_implementation() == 'PyPy':
    import pypysam
else:
    import pysam

import sortedcontainers
import networkx as nx
from configutator import ConfigMap
import io, itertools, collections
from sys import maxsize

from FamilyRecord import FamilyRecord

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

def familyName(record: pysam.AlignedSegment, barcode_mask: str) -> str:
    """
    Families are collapsed based on position, forward/reverse, +/-, barcode
    :param record:
    :param barcode_mask:
    :return:
    """
    family = "R" if record.is_reverse else "F"  # Keep forward and reverse reads separate
    bc = record.get_tag("BC")
    for i in range(len(barcode_mask)):
        if barcode_mask[i] == "1":
            family += bc[i]
    return family

class Families:
    def __init__(self, barcode_mask: str, barcode_distance: int):
        self._familyRecords = sortedcontainers.SortedDict()
        self._barcode_mask = barcode_mask
        self._barcode_distance = barcode_distance

    def addRecord(self, record: pysam.AlignedSegment):
        pos = self._familyRecords.get(record.reference_start)
        family = familyName(record, self._barcode_mask)
        if pos is None:
            self._familyRecords[record.reference_start] = {family: FamilyRecord(family, record.reference_start, record)}
        else:
            familyRecord = pos.get(family)
            if familyRecord:
                familyRecord.aggregate(record)
            else:
                pos[family] = FamilyRecord(family, record.reference_start, record)

    def _cutTreeToThreshold(self, forest, tree, thresholdDiameter):
        eccentricity = nx.eccentricity(tree, sp=forest['paths'])
        if thresholdDiameter <= nx.diameter(tree, eccentricity):
            peripheral = nx.periphery(tree, eccentricity)
            candidateEdges = {}
            for node1 in peripheral:
                """Determine all candidate edges to remove that will result in a peripheral node being within threshold"""
                for node2 in peripheral:
                    path = nx.bidirectional_dijkstra(tree, node1, node2)
                    length = 0
                    for node in path:
                        data = tree.get_edge_data(node1, node)
                        length += data['weight']
                        if length > thresholdDiameter:
                            candidateEdges[(node1, node)] = data
                            break
            maxRankEdge = None
            for edge, data in candidateEdges.items():
                """Determine candidate edge that will result in maximum total node weight"""
                tree.remove_edge(*edge)
                for component in nx.connected_components(tree):
                    rank = 0
                    for node in component:
                        rank += forest[node]['weight']
                    if rank > data['rank']:
                        data['rank'] = rank
                if maxRankEdge is None or maxRankEdge[1] < data['rank']:
                    maxRankEdge = (edge, data['rank'])
                tree.add_edge(*edge, data)
            forest.remove_edge(*maxRankEdge[0]) #Remove edge from forest to store global result
            tree.remove_edge(*maxRankEdge[0])
            for component in nx.connected_component_subgraphs(tree, False):
                """Recurse for each new subtree to ensure it is within threshold"""
                self._cutTreeToThreshold(forest, component, thresholdDiameter)

    def _collapsePosition(self, pos: dict):
        #Collapse to weighted minimum spanning forest with trees that have no path greater than 2*threshold maximising family weight. What to do with sequence validation?
        if pos:
            keys = pos.keys()
            graph = nx.empty_graph(len(pos)) #type: nx.Graph
            for node, data in graph.nodes_iter(True):
                data['weight'] = pos[keys[node]].size
            for u, v in itertools.combinations(range(len(pos)), 2):
                dist = hamming(keys[u], keys[v], self._barcode_distance)
                if dist < self._barcode_distance: graph.add_edge(u, v, weight=dist)
            thresholdDiameter = 2 * self._barcode_distance
            allShortestPaths = nx.all_pairs_dijkstra_path_length(graph)
            eccentricity = nx.eccentricity(graph, sp=allShortestPaths)
            if thresholdDiameter <= nx.diameter(graph, eccentricity):
                forest = nx.minimum_spanning_tree(graph) #type: nx.Graph
                forest['paths'] = allShortestPaths
                #TODO correct greedy algorithm for minimum eccentricity
                self._cutTreeToThreshold(forest, forest, thresholdDiameter)
            for component in nx.connected_components(forest):
                for node in component[1:]:
                    pos[keys[component[0]]] += pos[keys[node]]

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

    def __delitem__(self, pos):
        del self._familyRecords[pos]


@ConfigMap(inStream=None, outStream=None)
def collapse(inStream: io.IOBase, outStream: io.IOBase, barcode_distance: int, barcode_mask: str):
    """
    Collapse reads with the same start position, forward/reverse, and barcode into consensus record.
    :param inStream: A file or stream handle to read input data
    :param outStream: A file or stream handle to output data
    :param barcode_distance: The maximum number of differences allowed to be considered the same family
    :param barcode_mask: String mask to denote the positions within a barcode to include in the family name
    :return: None
    """
    inFile = pysam.AlignmentFile(inStream)
    outFile = pysam.AlignmentFile(outStream, "rbu")
    families = Families(barcode_mask, barcode_distance)
    startPos = 0
    for record in inFile.fetch(until_eof=True):
        if startPos != record.reference_start:
            for _, family in families.getFamiliesAt(startPos):
                outFile.write(family.toPysam())
            del families[startPos]
            startPos = record.reference_start
        families.addRecord(record)

if __name__ == "__main__":
    from sys import stdout, stdin, stderr, argv, maxsize, maxsize
    from configutator import loadConfig
    for argmap, paths in loadConfig(argv, (collapse,), title="Collapse V1.0"):
        collapse(paths[0] if len(paths) else stdin, open(paths[1], 'wb+') if len(paths) > 1 else stdout, **argmap[collapse])
        break