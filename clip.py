#!/usr/bin/env python3
import pysam
from sys import maxsize
from CigarIterator import CigarIterator, appendOrInc

default_cost = {
    pysam.CMATCH: lambda x: -x,
    pysam.CEQUAL: lambda x: -x,
    pysam.CDIFF: lambda x: -x,
    pysam.CREF_SKIP: lambda x: -x,
    pysam.CINS: lambda x: 6 + 1*(x-1),
    pysam.CDEL: lambda x: 3 + 1*(x-1)
}

def calculateAlignmentCost(record: pysam.AlignedSegment, start: int = 0, end: int = maxsize, costs = default_cost) -> int:
    cost = 0
    if start < record.reference_start:
        start = record.reference_start
    #if end > record.reference_end:
    #    end = record.reference_end
    if start > end:
        #No works needs to be done
        return cost
    i = CigarIterator(record)
    i.skipToRefPos(start)
    if not i.clipped: cost += costs[i.op](i.opEnd - i.opPos + 1)
    while not i.clipped and i.refPos+i.opLength-1 <= end:
        cost += costs[i.op](i.opLength)
        if not i.nextOp(): break
    if i.valid and not i.clipped:
        cost += costs[i.op](end - i.refPos + 1)
    return cost

def calculateMappingQuality(record: pysam.AlignedSegment) -> int:
    pass #TODO

def trimRecord(record: pysam.AlignedSegment, mate: pysam.AlignedSegment, start: int = 0, end: int = maxsize):
    if start > end:
        #No works needs to be done
        return
    if start < record.reference_start:
        start = record.reference_start
    #if end >= record.reference_end:
    #    end = record.reference_end-1
    ops = []
    nextMatch = None # Retain the first match after clipping to avoid setting the new start position to a deletion
    i = CigarIterator(record)
    hard = i.skipClipped(True)
    if hard: ops += [(pysam.CHARD_CLIP, hard)]
    if i.skipToRefPos(start):
        appendOrInc(ops, [pysam.CSOFT_CLIP, i.seqPos])
        appendOrInc(ops, [i.op, i.opEnd - i.opPos +1])
        if i.op in (pysam.CMATCH, pysam.CEQUAL, pysam.CDIFF): # Are we there yet?
            nextMatch = i.refPos
        dist = end - i.refPos
        while dist > i.opLength and i.nextOp():
            if nextMatch is None and i.op in (pysam.CMATCH, pysam.CEQUAL, pysam.CDIFF): # Are we there yet?
                nextMatch = i.refPos
            appendOrInc(ops, i.ops[i.opsI])
            if i.inRef: dist -= i.opLength
        if i.valid and dist <= i.opLength:
            if nextMatch is None and i.op in (pysam.CMATCH, pysam.CEQUAL, pysam.CDIFF): # Are we there yet?
                nextMatch = i.refPos
            appendOrInc(ops, [i.op, dist])
            appendOrInc(ops, [pysam.CSOFT_CLIP, i.opLength - dist])
        appendOrInc(ops, [pysam.CSOFT_CLIP, i.record.query_length - i.seqPos])
    else:
        #Soft clip entire read
        appendOrInc(ops, [pysam.CSOFT_CLIP, record.query_length])
        nextMatch = start
    record.cigartuples = ops
    record.reference_start = nextMatch
    # TODO update mapping quality
    #TODO rewrite MD
    mate.next_reference_start = record.reference_start

def mergeRecord(fromRecord: pysam.AlignedSegment, toRecord: pysam.AlignedSegment, refStart: int = -1, refEnd: int = maxsize, costs = default_cost) -> list:

    rStart = fromRecord.reference_start if fromRecord.reference_start > toRecord.reference_start else toRecord.reference_start
    rEnd = (fromRecord.reference_end if fromRecord.reference_end < toRecord.reference_end else toRecord.reference_end) - 1

    if refStart < rStart:
        refStart = rStart
    if refEnd > rEnd:
        refEnd = rEnd

    if fromRecord.reference_length == 0 or toRecord.reference_length == 0 or refStart >= refEnd:
        #No work needs to be done
        return

    ops = []
    toItr = CigarIterator(toRecord)

    #Copy in unaffected ops in non-overlapping region
    dist = refStart - toItr.refStart
    if not toItr.nextOp(): return
    while dist > toItr.opLength:
        appendOrInc(ops, list(toItr.opRange))
        if toItr.inRef: dist -= toItr.opLength
        if not toItr.nextOp(): return

    appendOrInc(ops, [toItr.op, dist])
    toItr.step(dist)
    seq = toRecord.query_sequence[:toItr.seqPos]
    qual = list(toRecord.query_qualities[:toItr.seqPos])

    fromItr = CigarIterator(fromRecord)
    if not fromItr.skipToRefPos(refStart): return

    toOptimal = None

    while toItr.refPos <= refEnd or fromItr.refPos <= refEnd:
        if toItr.op == fromItr.op:
            if toItr.inSeq:  #.getOp() in [pysam.CMATCH, pysam.CEQUAL, pysam.CDIFF]:
                if toItr.baseQual == fromItr.baseQual:
                    r = toItr if not toItr.matchesRef else fromItr
                    seq += r.seqBase #Keep the variant if possible
                    qual += [r.baseQual if r.matchesRef or toItr.seqBase == fromItr.seqBase else 3] # 3 = 50% probability of either base being correct
                elif toItr.baseQual > fromItr.baseQual:
                    seq += toItr.seqBase
                    qual += [toItr.baseQual]
                else:
                    seq += fromItr.seqBase
                    qual += [fromItr.baseQual]
            appendOrInc(ops, [toItr.op, 1])
        else:
            if toOptimal == None:
                # Dont calculate costs if unnecessary
                toOptimal = calculateAlignmentCost(toRecord, refStart, refEnd, costs) > calculateAlignmentCost(fromRecord, refStart, refEnd, costs)
            if toOptimal:
                if toItr.inSeq:
                    seq += toItr.seqBase
                    qual += [toItr.baseQual]
                appendOrInc(ops, [toItr.op, 1])
            else:
                if fromItr.inSeq:
                    seq += fromItr.seqBase
                    qual += [fromItr.baseQual]
                appendOrInc(ops, [fromItr.op, 1])
        if not fromItr.next():
            while toItr.next():  # Copy remainder of toRecord
                if toItr.inSeq:
                    seq += toItr.seqBase
                    qual += [toItr.baseQual]
                appendOrInc(ops, [toItr.op, 1])
            break
        if not toItr.next():
            break
    while toItr.next():  # Copy remainder of toRecord
        if toItr.inSeq:
            seq += toItr.seqBase
            qual += [toItr.baseQual]
        appendOrInc(ops, [toItr.op, 1])
    toRecord.cigartuples = ops
    toRecord.query_sequence = seq
    toRecord.query_qualities = qual
    # TODO update mapping quality

if __name__ == '__main__':
    import getopt
    from sys import stdout, stdin, stderr, argv
    threads = 1
    max_read_len = 1000
    outFormat = 'w'
    bufferingMode = 'p'

    #Command line options
    stderr.write("Clip Overlap v1.0\n")
    ops, paths = getopt.gnu_getopt(argv[1:], 't:m:ho:s:')
    for op, val in ops:
        if op == '-h':
            stderr.write("Clip overlapping reads from SAM/BAM/CRAM file\n"
                         "Use: clip.py [-t threads] [-m max_read_length] [input file path | < infile > outfile] [output file path]\n"
                         "If no paths are given stdin and stdout are used.\n"
                         "-t # of threads to use for processing (Default=1)\n"
                         "-m Maximum possible read length in data (Default=1000)\n"
                         "-o [sbuc] Output format: s=SAM (Default), b=BAM compressed, bu=BAM uncompressed, c=CRAM\n"
                         "-s [mcp] Buffering mode:\n"
                         "\tm=Maintain input order\n"
                         "\tc=Maintain coordinate order (Requires sorted input)\n"
                         "\tp=Output in arbitrary order (Minimal RAM, Default)\n")
            exit()
        elif op == '-t':
            threads = int(val or 1)
        elif op == '-m':
            max_read_len = int(val or 1000)
        elif op == '-o' and val != 's':
            outFormat += val
        elif op == '-s':
            bufferingMode = val
    else:
        stderr.write("Use -h for help.\n")

    import multiprocessing

    #from collections import deque #TODO comment out for multithreading
    pool = []
    #class qwrap(deque): #TODO comment out for multithreading
    #    def get(self):
    #        return self.popleft()
    #    def put(self, item):
    #        return self.append(item)
    #    def empty(self):
    #        return len(self) == 0
    #    def task_done(self):
    #        pass

    queueIn = multiprocessing.JoinableQueue(3*threads)
    queueOut = multiprocessing.JoinableQueue(3*threads)
    #queueIn = qwrap()
    #queueOut = qwrap()
    mateBuffer = {}
    inFile = pysam.AlignmentFile(paths[0] if len(paths) else stdin)
    #outFile = pysam.AlignmentFile(paths[1] if len(paths) > 1 else stdout, outFormat, template=inFile) #TODO comment out for threading
    def worker():
        while True: #TODO Uncomment for threading
            if not queueIn.empty():
                firstRecord, secondRecord = queueIn.get()
                if firstRecord[1].reference_start < secondRecord[1].reference_start:
                    leftRecord = firstRecord[1]
                    rightRecord = secondRecord[1]
                else:
                    rightRecord = firstRecord[1]
                    leftRecord = secondRecord[1]
                mergeRecord(leftRecord, rightRecord)
                trimRecord(leftRecord, rightRecord, rightRecord.reference_end)
                queueOut.put(firstRecord)
                queueOut.put(secondRecord)
                queueIn.task_done()

    def writer():
        outFile = pysam.AlignmentFile(paths[1] if len(paths) > 1 else stdout, outFormat, template=inFile)
        if bufferingMode == 'm':
            from sortedcontainers import SortedListWithKey
            writeBuffer = SortedListWithKey(key=lambda x: x[0])
            nextIndex = 0 #type: int
            while True:
                record = queueOut.get()
                if record[0] == nextIndex:
                    nextIndex += 1
                    outFile.write(record[1])
                    queueOut.task_done()
                else:
                    writeBuffer.add(record)
                if len(writeBuffer) and writeBuffer[0][0] == nextIndex:
                    record = writeBuffer.pop(0)
                    nextIndex += 1
                    outFile.write(record[1])
                    queueOut.task_done()
        elif bufferingMode == 'c':
            from sortedcontainers import SortedListWithKey
            sortBuffer = SortedListWithKey(key=lambda x: x[0]) #type: List[(int, pysam.AlignedSegment, int)]
            writeBuffer = SortedListWithKey(key=lambda x: x[1].reference_start) #type: List[(int, pysam.AlignedSegment, int)]
            nextIndex = 0  # type: int
            while True:
                sortBuffer.add(queueOut.get())
                if len(sortBuffer) and sortBuffer[0][0] == nextIndex:
                    writeBuffer.add(sortBuffer.pop(0))
                    nextIndex += 1
                if len(writeBuffer) > 2 and writeBuffer[-1][2] > writeBuffer[0][1].reference_start:
                    record = writeBuffer.pop(0)
                    outFile.write(record[1])
                    queueOut.task_done()
        else:
            while True:
                record = queueOut.get()
                outFile.write(record[1])
                queueOut.task_done()

    pool.append(multiprocessing.Process(target=writer))  #TODO uncomment for threading
    for _ in range(threads):
        pool.append(multiprocessing.Process(target=worker))
    for p in pool:
        p.start()

    #TODO Comment out for threading
    #from sortedcontainers import SortedListWithKey
    #if bufferingMode == 'm':
    #    writeBuffer = SortedListWithKey(key=lambda x: x[0])
    #elif bufferingMode == 'c':
    #    writeBuffer = SortedListWithKey(key=lambda x: x[1].reference_start)  # type: List[(int, pysam.AlignedSegment)]
    #
    #nextIndex = 0  # type: int
    #---

    i = 0 # Tracks the order the records were read
    for record in inFile:
        #Skip clipping code if the records don't overlap
        #TODO Uncomment for threading
        if record.is_paired and record.reference_start-max_read_len < record.next_reference_start < record.reference_end:
            secondRecord = mateBuffer.get(record.query_name, None)
            if not secondRecord:
                mateBuffer[record.query_name] = (i, record, record.reference_start)
            else:
                queueIn.put(((i, record, record.reference_start), secondRecord))
        else:
            queueOut.put((i, record, record.reference_start))
        #---
        #TODO Comment out for treading
        #if record.is_paired and record.reference_start - max_read_len < record.next_reference_start < record.reference_end:
        #    firstRecord = mateBuffer.get(record.query_name, None)
        #    if not firstRecord:
        #        mateBuffer[record.query_name] = (i if bufferingMode == 'm' else record.reference_start, record)
        #    else:
        #        queueIn.put((firstRecord, (i if bufferingMode == 'm' else record.reference_start, record)))
        #        worker()
        #        if bufferingMode == 'p':
        #            outFile.write(queueOut.get()[1])
        #            outFile.write(queueOut.get()[1])
        #        else:
        #            writeBuffer.add(queueOut.get())
        #            writeBuffer.add(queueOut.get())
        #else:
        #    writeBuffer.add((i if bufferingMode == 'm' else record.reference_start, record))
        #if bufferingMode == 'm':
        #    while len(writeBuffer) and writeBuffer[0][0] == nextIndex:
        #        r = writeBuffer.pop(0)
        #        nextIndex += 1
        #        outFile.write(r[1])
        #elif bufferingMode == 'c':
        #    while len(writeBuffer) > 2 and writeBuffer[-1][0] > writeBuffer[0][1].reference_start and writeBuffer[-1][0] > writeBuffer[0][0]:
        #        r = writeBuffer.pop(0)
        #        outFile.write(r[1])
        #---
        i += 1

    queueIn.join() #TODO uncomment for threading
    #while not queueIn.empty():#TODO comment out for threading
    #    worker()
    for _, record in mateBuffer.items():
        queueOut.put(record)
    queueOut.join() #TODO uncomment for threading
    #if bufferingMode != 's':
    #    for record in writeBuffer:#TODO comment out for threading
    #        outFile.write(record[1])#TODO comment out for threading
    inFile.close()
    #outFile.close()