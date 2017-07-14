#!/usr/bin/env python3
import pysam
from sys import maxsize
from CigarIterator import CigarIterator, appendOrInc
import multiprocessing, ctypes
import parapysam

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
    # TODO rewrite MD
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

class WorkerProcess(multiprocessing.Process):
    def __init__(self, header=None):
        super().__init__()
        self.start(header)

    def start(self, header):
        self.orderPipeR, self.orderPipeW = multiprocessing.Pipe(False)
        self.recordsInR, self.recordsInW = parapysam.RecordPipe(header)
        self.recordsOutR, self.recordsOutW = parapysam.RecordPipe(header)
        super().start()

    def run(self):
        try:
            firstRecord = next(self.recordsInR)
            secondRecord = next(self.recordsInR)
            while firstRecord and secondRecord:
                if firstRecord.reference_start < secondRecord.reference_start:
                    leftRecord = firstRecord
                    rightRecord = secondRecord
                else:
                    rightRecord = firstRecord
                    leftRecord = secondRecord
                mergeRecord(leftRecord, rightRecord)
                trimRecord(leftRecord, rightRecord, rightRecord.reference_end)
                self.recordsOutW.write(firstRecord)
                self.recordsOutW.write(secondRecord)
                firstRecord = next(self.recordsInR)
                secondRecord = next(self.recordsInR)
        except StopIteration:
            pass
        self.recordsOutW.close()

    def stop(self):
        self.recordsInW.close()
        self.join()

    def __del__(self):
        self.stop()
        #super().__del__()

    def sendMatePair(self, index1, record1, index2, record2):
        self.recordsInW.write(record1)
        self.orderPipeW.send(index1)
        self.recordsInW.write(record2)
        self.orderPipeW.send(index2)

class WriterProcess(multiprocessing.Process):
    def __init__(self, outFH, outFormat, pool, header, ordered):
        super().__init__()
        self.pool = pool
        self.ordered = ordered
        self.outFH = outFH
        self.outFormat = outFormat
        self.nextIndex = multiprocessing.RawValue(ctypes.c_long, 0)
        self.bufferedCount = multiprocessing.RawValue(ctypes.c_long, 0)
        self.start(header)

    def start(self, header):
        self.header = header
        self.orderPipeR, self.orderPipeW = multiprocessing.Pipe(False)
        self.recordsOutR, self.recordsOutW = parapysam.RecordPipe(header)
        self.recordsOutW.flush()
        super().start()

    def run(self):
        outFile = pysam.AlignmentFile(self.outFH, self.outFormat, header=self.header)
        if self.ordered:
            self.writeOrdered(outFile)
        else:
            self.write(outFile)
        outFile.close()

    def stop(self):
        self.recordsOutW.close()
        for p in self.pool:
            p.stop()
        self.join()
        self.orderPipeR.close()

    def __del__(self):
        self.stop()
        #super().__del__()

    def sendUnpaired(self, index, record):
        self.recordsOutW.write(record)
        self.orderPipeW.send(index)

    def writeOrdered(self, outFile):
        from sortedcontainers import SortedListWithKey
        writeBuffer = SortedListWithKey(key=lambda x: x[0])
        self.nextIndex.value = 0  # type: int
        running = True
        while running:
            running = False
            for p in self.pool + [self]: #type: WorkerProcess
                if p.recordsOutR.eof:
                    continue
                running = True
                if not p.recordsOutR.poll():
                    continue
                result = (p.orderPipeR.recv(), p.recordsOutR.read())
                if result[0] == self.nextIndex.value:
                    self.nextIndex.value += 1
                    outFile.write(result[1])
                else:
                    writeBuffer.add(result)
                    self.bufferedCount.value += 1
                while len(writeBuffer) and writeBuffer[0][0] == self.nextIndex.value:
                    result = writeBuffer.pop(0)
                    self.bufferedCount.value -= 1
                    self.nextIndex.value += 1
                    outFile.write(result[1])

    def write(self, outFile):
        running = True
        while running:
            running = False
            for p in self.pool + [self]: #type: WorkerProcess
                if p.recordsOutR.closed:
                    continue
                running = True
                if not p.recordsOutR.poll():
                    continue
                self.nextIndex.value = p.orderPipeR.recv()
                outFile.write(p.recordsOutW.read())

def status(mateCount, bufferedCount, nextIndex):
    import time
    status = ''
    while True:
        status = "\rMates buffered: {:>10}\tPending output: {:>10}\tWaiting for read number: {:>10}\t".format(mateCount.value, bufferedCount.value, nextIndex.value)
        stderr.write(status)
        time.sleep(1)

def clip(paths, threads, maxTLen, outFormat, ordered, verbose):
    pool = []
    mateBuffer = {}
    inFile = pysam.AlignmentFile(paths[0] if len(paths) else stdin)
    inFileItr = inFile.fetch(until_eof=True)

    # Begin processing
    for _ in range(threads):
        pool.append(WorkerProcess(inFile.header))

    writer = WriterProcess(open(paths[1], 'wb+') if len(paths) > 1 else stdout, outFormat, pool, inFile.header, ordered)

    def poolLooper():
        while True:
            for p in pool:
                yield p

    poolLoop = poolLooper()
    mateCount = multiprocessing.RawValue(ctypes.c_long, 0)

    statusProc = multiprocessing.Process(target=status, args=(mateCount, writer.bufferedCount, writer.nextIndex))
    if verbose:
        statusProc.start()

    i = 0  # Tracks the order the records were read
    for record in inFileItr:
        # Skip clipping code if the records don't overlap
        if not record.is_unmapped and not record.mate_is_unmapped and not record.is_supplementary and record.is_paired and abs(record.template_length) < maxTLen and record.reference_name == record.next_reference_name: #record.reference_start - maxTLen < record.next_reference_start < record.reference_end:
            secondRecord = mateBuffer.get(record.query_name, None)
            if not secondRecord:
                mateBuffer[record.query_name] = (i, record)
                mateCount.value += 1
            else:
                next(poolLoop).sendMatePair(i, record, *secondRecord)
                del mateBuffer[record.query_name]
                mateCount.value -= 1
        else:
            writer.sendUnpaired(i, record)
        i += 1
    inFile.close()

    for _, record in mateBuffer.items():
        mateCount.value -= 1
        writer.sendUnpaired(*record)
    writer.stop()
    if verbose:
        statusProc.terminate()

if __name__ == '__main__':
    import getopt
    from sys import stdout, stdin, stderr, argv
    threads = 1
    maxTLen = 1000
    outFormat = 'w'
    ordered = False
    verbose = False

    #Command line options
    stderr.write("Clip Overlap v1.0\n")
    ops, paths = getopt.gnu_getopt(argv[1:], 't:m:ho:sv')
    if not len(ops) and not len(paths):
        stderr.write("Use -h for help.\n")
    else:
        for op, val in ops:
            if op == '-h':
                stderr.write("Clip overlapping reads from SAM/BAM/CRAM file\n"
                             "Use: clip.py [-tmos] [input file path | < infile > outfile] [output file path]\n"
                             "If no paths are given stdin and stdout are used.\n"
                             "-t # Threads to use for processing (Default=1)\n"
                             "-m # Maximum template length guaranteeing no read overlap (Default=1000)\n"
                             "-o [sbuc] Output format: s=SAM (Default), b=BAM compressed, bu=BAM uncompressed, c=CRAM\n"
                             "-s Maintain input order (High depth regions may fill RAM), if not set will output in arbitrary order (Minimal RAM)\n"
                             "-v Verbose status output\n")
                exit()
            elif op == '-t':
                threads = int(val or 1)
            elif op == '-m':
                maxTLen = int(val or 1000)
            elif op == '-o' and val != 's':
                outFormat += val
            elif op == '-s':
                ordered = True
            elif op == '-v':
                verbose = True

    clip(paths, threads, maxTLen, outFormat, ordered, verbose)