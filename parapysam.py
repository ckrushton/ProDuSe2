import multiprocessing, pysam, threading

class RecordPipe:
    _flushRecordName = "___%!%DummyRecordShouldNotBeInOutput%!%___"
    def __new__(cls, header):
        r, w = multiprocessing.Pipe(False)
        closeEvent = multiprocessing.Event()
        return RecordOutPipe(r, closeEvent), RecordInPipe(w, closeEvent, header)

class RecordInPipe:
    # This is to deal with the fact that pysam locks up when waiting for input but gets EOF
    _flushRecord = pysam.AlignedSegment()
    _flushRecord.query_name = RecordPipe._flushRecordName
    def __init__(self, pipe, closeEvent, header):
        self._pipe = pipe #type: multiprocessing.Connection
        self._closeEvent = closeEvent
        self._header = header
        self._hts = None

    def close(self):
        if self._closeEvent.is_set():
            return
        self.flush()
        if self._hts:
            self._hts.close()
        self._pipe.close()
        self._closeEvent.set()

    @property
    def closed(self):
        return self._closeEvent.is_set()

    def __del__(self):
        self.close()

    def flush(self):
        if not self._hts:
            self._hts = pysam.AlignmentFile(self._pipe, 'wbu', header=self._header)
        for i in range(30):
            self._hts.write(self._flushRecord)

    def write(self, record):
        if not self._hts:
            self._hts = pysam.AlignmentFile(self._pipe, 'wbu', header=self._header)
            self.flush()
        self._hts.write(record)

class RecordOutPipe:
    def __init__(self, pipe, closeEvent):
        self._pipe = pipe
        self._closeEvent = closeEvent
        self._hts = None
        self._htsItr = None
        self._nextRecord = None
        self._thread = None
        self._recordAvailable = threading.Event()
        self._getNextRecord = threading.Event()

    def _initHTS(self):
        if not self._thread:
            self._thread = threading.Thread(target=self._read)
            self._thread.start()

    def close(self):
        if self._hts:
            self._hts.close()
        self._pipe.close()

    @property
    def closed(self):
        return self._closeEvent.is_set()

    @property
    def eof(self):
        return self.closed and not self.poll() and self._nextRecord and self._nextRecord.query_name == RecordPipe._flushRecordName

    def __del__(self):
        self.close()

    def _read(self):
        if not self._htsItr:
            self._hts = pysam.AlignmentFile(self._pipe.fileno(), check_header=False, check_sq=False)
            self._htsItr = self._hts.fetch(until_eof=True)
        while not self.eof:
            self._nextRecord = next(self._htsItr)
            while self._nextRecord.query_name == RecordPipe._flushRecordName:
                self._nextRecord = next(self._htsItr)
            self._recordAvailable.set()
            self._getNextRecord.wait()
            self._getNextRecord.clear()

    def poll(self):
        self._initHTS()
        return self._recordAvailable.is_set()

    def read(self):
        self._initHTS()
        self._getNextRecord.set()
        self._recordAvailable.wait()
        self._recordAvailable.clear()
        return self._nextRecord

    def __iter__(self):
        return self

    def __next__(self):
        self._initHTS()
        self._getNextRecord.set()
        while True:
            if self._recordAvailable.wait(0.1):
                self._recordAvailable.clear()
                return self._nextRecord
            elif self.closed:
                self.close()
                raise StopIteration