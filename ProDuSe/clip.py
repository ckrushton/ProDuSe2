#!/usr/bin/env python3
import platform
if platform.python_implementation() == 'PyPy':
    import pypysam
else:
    import pysam
from sys import maxsize, stderr
import multiprocessing, ctypes, io

# If running as a script, this works fine

try:
    from CigarIterator import CigarIterator, appendOrInc
    import parapysam
# If installed
except ModuleNotFoundError:
    from ProDuSe.CigarIterator import CigarIterator, appendOrInc
    from ProDuSe import parapysam

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
        #No work needs to be done
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
    # Mapping quality isn't updated as it reflects the quality of the original alignment
    return record.mapping_quality

def trimRecord(record: pysam.AlignedSegment, mate: pysam.AlignedSegment, start: int = 0, end: int = maxsize):
    if start > end:
        #No works needs to be done
        return
    if start < record.reference_start:
        start = record.reference_start
    #if end >= record.reference_end:
    #    end = record.reference_end-1
    ops = []

    # Retain the first match after clipping to avoid setting the new start position to a deletion
    nextMatch = None
    i = CigarIterator(record)
    def checkForFirstMatch():
        nonlocal nextMatch, i
        if nextMatch is None and i.op in (pysam.CMATCH, pysam.CEQUAL, pysam.CDIFF):
            nextMatch = i.refPos

    #Jump past hard clipping and retain op
    hard = i.skipClipped(True)
    if hard: ops += [(pysam.CHARD_CLIP, hard)]
    #Skip to reference starting position
    if start < end and i.skipToRefPos(start):
        #Soft clip everything before start
        appendOrInc(ops, [pysam.CSOFT_CLIP, i.seqPos])
        dist = end - i.refPos #Calculate reference distance remaining to end

        #Copy in all ops before end
        while dist > i.opRemaining or not i.inRef:
            checkForFirstMatch()
            appendOrInc(ops, [i.op, i.opRemaining])
            if i.inRef: dist -= i.opRemaining
            if not i.nextOp(): break

        #If end within op, copy in remainder
        if i.valid:
            checkForFirstMatch()
            appendOrInc(ops, [i.op, dist])
            if not i.inSeq: dist = 0
            appendOrInc(ops, [pysam.CSOFT_CLIP, i.record.query_length - i.seqPos - dist])
            # Retain hard clip at end
            if len(i.ops) and i.ops[-1][0] == pysam.CHARD_CLIP:
                appendOrInc(ops, i.ops[-1])
    else:
        #Soft clip entire read
        appendOrInc(ops, [pysam.CSOFT_CLIP, record.query_length])
        #Retain hard clip at end
        if len(i.ops) and i.ops[-1][0] == pysam.CHARD_CLIP:
            appendOrInc(ops, i.ops[-1])

    #Update record
    record.cigartuples = ops
    if nextMatch is not None:
        record.set_tag('OS', record.reference_start)
        record.reference_start = nextMatch
    record.mapping_quality = calculateMappingQuality(record)
    # TODO rewrite MD
    mate.next_reference_start = record.reference_start

def mergeRecord(fromRecord: pysam.AlignedSegment, toRecord: pysam.AlignedSegment, refStart: int = -1, refEnd: int = maxsize, costs = default_cost) -> list:
    overlapStart = fromRecord.reference_start if fromRecord.reference_start > toRecord.reference_start else toRecord.reference_start
    overlapEnd = (fromRecord.reference_end if fromRecord.reference_end < toRecord.reference_end else toRecord.reference_end) - 1

    if refStart < overlapStart:
        refStart = overlapStart
    if refEnd > overlapEnd:
        refEnd = overlapEnd

    if fromRecord.reference_length == 0 or toRecord.reference_length == 0 or refStart >= refEnd:
        #No work needs to be done
        return

    ops = []
    toItr = CigarIterator(toRecord)

    #Copy in unaffected ops in non-overlapping region
    dist = refStart - toItr.refStart
    if not toItr.nextOp(): return
    while not toItr.inRef or dist > toItr.opLength:
        appendOrInc(ops, list(toItr.opRange))
        if toItr.inRef: dist -= toItr.opLength
        if not toItr.nextOp(): return

    appendOrInc(ops, [toItr.op, dist])
    toItr.step(dist)
    seq = toRecord.query_sequence[:toItr.seqPos]
    qual = list(toRecord.query_qualities[:toItr.seqPos])

    #Init fromRecord iterator
    fromItr = CigarIterator(fromRecord)
    if not fromItr.skipToRefPos(refStart): return

    toOptimal = None

    while toItr.refPos <= refEnd and fromItr.refPos <= refEnd:
        if toItr.op == fromItr.op:
            if toItr.inSeq:
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
    while toItr.valid:  # Copy remainder of toRecord
        if toItr.inSeq:
            seq += toItr.seqBase
            qual += [toItr.baseQual]
        appendOrInc(ops, [toItr.op, 1])
        toItr.next()
    toRecord.cigartuples = ops
    toRecord.query_sequence = seq
    toRecord.query_qualities = qual
    # TODO update mapping quality
    # TODO Store tag for which strand the base originated from

class WorkerProcess(parapysam.OrderedWorker):
    def __init__(self, alternate, clipOnly):
        self.alternate = alternate
        self.clipOnly = clipOnly
        super().__init__()

    def work(self):
        flipflop = False
        try:
            firstRecord = self.receiveRecord()
            secondRecord = self.receiveRecord()
            while firstRecord and secondRecord:
                if firstRecord.reference_start < secondRecord.reference_start:
                    leftRecord = firstRecord
                    rightRecord = secondRecord
                else:
                    rightRecord = firstRecord
                    leftRecord = secondRecord
                # alternate clipped read to avoid strand bias
                if self.alternate and flipflop:
                    if not self.clipOnly: mergeRecord(rightRecord, leftRecord)
                    trimRecord(rightRecord, leftRecord, leftRecord.reference_end-1)
                    flipflop = False
                else:
                    if not self.clipOnly: mergeRecord(leftRecord, rightRecord)
                    trimRecord(leftRecord, rightRecord, 0, rightRecord.reference_start)
                    flipflop = True
                self.sendRecord(firstRecord)
                self.sendRecord(secondRecord)
                firstRecord = self.receiveRecord()
                secondRecord = self.receiveRecord()
        except StopIteration:
            pass

class WriterProcess(parapysam.OrderedWorker):
    def __init__(self, outFH, outFormat, pool, ordered, maxTLen):
        super().__init__()
        self.pool = pool
        self.ordered = ordered
        self.outFH = outFH
        self.outFormat = outFormat
        self.nextIndex = multiprocessing.RawValue(ctypes.c_long, 0)
        self.bufferedCount = multiprocessing.RawValue(ctypes.c_long, 0)
        self.maxTLen = maxTLen

    def work(self):
        outFile = pysam.AlignmentFile(self.outFH, 'w' + self.outFormat, header=self.header)
        if self.ordered:
            self.writeOrdered(outFile)
        else:
            self.write(outFile)
        outFile.close()

    def stop(self):
        super().stop(False)
        for p in self.pool:
            p.stop()
        self.join()

    def writeOrdered(self, outFile):
        import heapq
        class HeapNode:
            __slots__ = 'key', 'value'
            def __init__(self, key, value):
                self.key, self.value = key, value
            def __lt__(self, other):
                return self.key < other.key
        writeBuffer = []
        self.nextIndex.value = 0  # type: int
        running = True
        while running:
            running = False
            for p in self.pool + [self]: #type: WorkerProcess
                if p.checkEOF():
                    continue
                running = True
                if not p.pollRecord():
                    continue
                result = p.receiveOrderedRecord()
                if self.ordered != 'c' and result[0] == self.nextIndex.value:
                    self.nextIndex.value += 1
                    outFile.write(result[1])
                else:
                    try:
                        startPos = result[1].get_tag('OS') #type: int
                    except KeyError:
                        startPos = result[1].reference_start
                    if self.ordered == 'c':
                        if self.nextIndex.value < startPos:
                            self.nextIndex.value = startPos
                        heapq.heappush(writeBuffer, HeapNode(startPos, result))
                    else:
                        heapq.heappush(writeBuffer, HeapNode(result[0], result))
                    self.bufferedCount.value += 1
                while len(writeBuffer) and (writeBuffer[0].key + self.maxTLen < self.nextIndex.value) if self.ordered == 'c' else (writeBuffer[0].key == self.nextIndex.value):
                    result = heapq.heappop(writeBuffer)
                    self.bufferedCount.value -= 1
                    if self.ordered != 'c': self.nextIndex.value += 1
                    outFile.write(result.value[1])

    def write(self, outFile):
        running = True
        while running:
            running = False
            for p in self.pool + [self]: #type: WorkerProcess
                if p.checkEOF():
                    continue
                running = True
                if not p.pollRecord():
                    continue
                self.nextIndex.value, record = p.receiveOrderedRecord()
                outFile.write(record)

def status(mateCount, bufferedCount, nextIndex, logStream: io.IOBase = stderr):
    import time
    logStream.write("\n") # This will be deleted by the next write
    while True:
        status = "\x1b[F\x1b[2K\rMates buffered: {:>10}\tPending output: {:>10}\tWaiting for read number: {:>10}\n".format(mateCount.value, bufferedCount.value, nextIndex.value)
        logStream.write(status)
        time.sleep(0.2)

def clip(inStream: io.IOBase, outStream: io.IOBase, threads: int = 8, maxTLen: int = 1000, outFormat:str = 'bu', ordered=False, alternate=False, trimBarcode=False, clipOnly=False, verbose=False, logStream: io.IOBase=stderr):
    pool = []
    mateBuffer = {}
    inFile = pysam.AlignmentFile(inStream)
    inFileItr = inFile.fetch(until_eof=True)

    # Begin processing
    for _ in range(threads):
        worker = WorkerProcess(alternate, clipOnly)
        worker.start(inFile.header)
        pool.append(worker)

    writer = WriterProcess(outStream, outFormat, pool, 'c' if alternate and ordered else ordered, maxTLen)
    writer.start(inFile.header)

    def poolLooper():
        while True:
            for p in pool:
                yield p

    poolLoop = poolLooper()
    mateCount = multiprocessing.RawValue(ctypes.c_long, 0)

    statusProc = multiprocessing.Process(target=status, args=(mateCount, writer.bufferedCount, writer.nextIndex, logStream))
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
            writer.sendOrderedRecord(i, record)
        i += 1
    inFile.close()

    for _, record in mateBuffer.items():
        mateCount.value -= 1
        writer.sendOrderedRecord(*record)
    writer.stop()
    if verbose:
        statusProc.terminate()
        logStream.write("\x1b[F\x1b[2K\rCompleted.\n")
    outStream.close()


def main(args=None):
    import getopt, os, errno
    from sys import stdout, stdin, argv
    threads = 2
    maxTLen = 1000
    outFormat = ''
    ordered = False
    alternate = False
    clipOnly = False
    trimBarcode = False
    verbose = False

    if args is None:
        args = argv[1:]
    #Command line options
    stderr.write("Clip v1.0\n")
    ops, paths = getopt.gnu_getopt(args, 't:m:ho:sva')
    if not len(ops) and not len(paths):
        stderr.write("Waiting on stdin. Try -h for help. Use ctrl+c to cancel.\n")

    else:
        for op, val in ops:
            if op == '-h':
                stderr.write("Clip overlapping reads from SAM/BAM/CRAM file\n"
                             "Use: clip.py [-tmos] [input file path | < infile > outfile] [output file path]\n"
                             "If no paths are given stdin and stdout are used.\n"
                             "-t # Threads to use for processing (Default=1)\n"
                             "-m # Maximum template length guaranteeing no read overlap (Default=1000)\n"
                             "-a Alternate strand being clipped to avoid strand bias (RAM intensive)\n"
                             "-b Trim trailing barcode region using BC tag"
                             "-c Clip only, do not merge clipped region into mate\n"
                             "-o [sbuc] Output format: s=SAM (Default), b=BAM compressed, bu=BAM uncompressed, c=CRAM\n"
                             "-s Maintain input order (High depth regions may fill RAM), if not set will output in arbitrary order (Minimal RAM)\n"
                             "-v Verbose status output\n")
                exit(1)
            elif op == '-t':
                threads = int(val or 1)
            elif op == '-m':
                maxTLen = int(val or 1000)
            elif op == '-a':
                alternate = True
            elif op == '-c':
                clipOnly = True
            elif op == '-b':
                trimBarcode = True
            elif op == '-o' and val != 's':
                outFormat += val
            elif op == '-s':
                ordered = True
            elif op == '-v':
                verbose = True

    if len(paths) > 1 and not os.path.exists(os.path.dirname(paths[1])):
        try:
            os.makedirs(os.path.dirname(paths[1]))
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

    clip(paths[0] if len(paths) else stdin, open(paths[1], 'wb+') if len(paths) > 1 else stdout, threads, maxTLen, outFormat, ordered, alternate, trimBarcode, verbose=verbose)

if __name__ == '__main__':
    main()
