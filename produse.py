#!/usr/bin/env python3
import sys, os, tempfile, shutil, gzip

from collapse import collapse
from clip import clip
from trim import trim

from multiprocessing import Process
from subprocess import Popen, PIPE

from configutator import ConfigMap, loadConfig
from valitator import Path, PathOrNone, Executable

class LoggedProcess(Process):
    def start(self, pipePath):
        tmp = sys.stderr
        sys.stderr = open(pipePath, 'wb')
        super().start()
        sys.stderr = tmp

@ConfigMap(fastq1='fastqs[0]', fastq2='fastqs[1]', output='join(``, [output, name, `.bam`])')
def ProDuSe(fastq1: Path, fastq2: Path, reference: Path, output: str, bwa: Executable, samtools: Executable, bed: PathOrNone = None, verbose: bool=False):
    """
    Manages the subprocesses and interprocess communication of Trim, BWA, Clip, samtools, Collapse, and Filter
    :param fastq1: Path to the first fastq
    :param fastq2: Path to the second fastq
    :param reference: Path to the reference fasta
    :param bed: Path to a bed file containing regions to restrict reads to (optional)
    :param output: Path to output file
    :param bwa: Path to bwa executable and any added parameters ex. mem or -t
    :param samtools: Path to samtools executable
    :return: None
    """
    global config
    stderr = sys.stderr
    temp_dir = tempfile.mkdtemp()
    stderr.write("Pipes created in {}\n".format(temp_dir))
    trimOut1Path = os.path.join(temp_dir, 'trimOut1')
    os.mkfifo(trimOut1Path)
    trimOut2Path = os.path.join(temp_dir, 'trimOut2')
    os.mkfifo(trimOut2Path)

    # Check if fastqs gzipped
    inFile1 = open(fastq1, 'rb')
    if inFile1.peek(2)[:2] == b'\037\213':
        inFile1 = gzip.open(inFile1)
    inFile2 = open(fastq2, 'rb')
    if inFile2.peek(2)[:2] == b'\037\213':
        inFile2 = gzip.open(inFile2)

    # Start bwa subprocess
    stderr.write("Starting BWA subprocess..\n")
    bwa = Popen(bwa.split(' ') + ['-C', reference, trimOut1Path, trimOut2Path], stdout=PIPE, stderr=PIPE)
    clipIn = bwa.stdout

    # Start trim subprocesses
    stderr.write("Starting Trim subprocesses..\n")
    trimDebug1 = os.path.join(temp_dir, 'trimStdErr1')
    trimDebug2 = os.path.join(temp_dir, 'trimStdErr2')
    os.mkfifo(trimDebug1)
    os.mkfifo(trimDebug2)
    trimDebug1In = open(trimDebug1, 'rb')
    trimDebug2In = open(trimDebug2, 'rb')
    trimProc1 = LoggedProcess(target=trim, args=(inFile1, open(trimOut1Path, 'wb')), kwargs=dict(verbose=verbose, **config[trim]))
    trimProc1.start(trimDebug1)
    trimProc2 = LoggedProcess(target=trim, args=(inFile2, open(trimOut2Path, 'wb')), kwargs=dict(verbose=verbose, **config[trim]))
    trimProc2.start(trimDebug2)

    view = None
    if bed:
        #Start samtools view to restrict to coordinates
        stderr.write("Starting samtools view subprocess..\n")
        view = Popen([samtools, "view", "-uL", bed, "-"], stdin=clipIn, stdout=PIPE, stderr=PIPE)
        clipIn = view.stdout

    # Start sort by coord
    stderr.write("Starting sort subprocess..\n")
    sort = Popen([samtools, "sort", "-l0", "-"], stdout=PIPE, stdin=PIPE, stderr=PIPE)

    #Start clipping subprocess
    stderr.write("Starting clipping subprocess..\n")
    clipDebug = os.path.join(temp_dir, 'clipStdErr')
    os.mkfifo(clipDebug)
    clipDebugIn = open(clipDebug, 'rb')
    clipProc = LoggedProcess(target=clip, args=(clipIn, sort.stdin), kwargs=dict(alternate=True, verbose=verbose))
    clipProc.start(clipDebug)

    stderr.write("Starting collapse subprocess..\n")
    collapseDebug = os.path.join(temp_dir, 'collapseStdErr')
    collapseDebugIn = open(collapseDebug, 'rb')
    collapseProc = LoggedProcess(target=collapse, kwargs=dict(inStream=sort.stdout, outStream=open(output, 'wb+'), verbose=verbose, **config[collapse]))
    collapseProc.start(collapseDebug)

    stderr.write("Writing to {}\n".format(output))
    trimBuffer1 = ''
    trimBuffer2 = ''
    bwaBuffer = []
    clipBuffer = ''
    collapseBuffer = ''
    while verbose and (trimProc1.is_alive() or trimProc2.is_alive() or clipProc.is_alive() or collapseProc.is_alive() or bwa.poll() or sort.poll()):
        stderr.write("\x1b[2J\x1b[H")  # Clear screen
        if trimDebug1In.peek(): trimBuffer1 = trimDebug1In.readline()
        if trimDebug2In.peek(): trimBuffer2 = trimDebug2In.readline()
        for _ in range(5):
            if not bwa.stderr.peek(): break
            bwaBuffer.pop(0)
            bwaBuffer.append(bwa.stderr.readline())
        if clipDebugIn.peek(): clipBuffer = clipDebugIn.readline()
        if collapseDebugIn.peek(): collapseBuffer = collapseDebugIn.readline()
        stderr.write("Writing to {}\nTrim:\n{}\n{}\n\nBWA:\n{}\n\nClip:\n{}\n\nCollapse:\n{}".format(output, trimBuffer1, trimBuffer2, ''.join(bwaBuffer), clipBuffer, collapseBuffer))

    trimProc1.join()
    trimProc2.join()
    clipProc.join()
    collapseProc.join()
    bwa.wait()
    if view:
        view.wait()
    sort.wait()
    inFile1.close()
    inFile2.close()
    shutil.rmtree(temp_dir)

if __name__ == "__main__":
    ConfigMap(_func='trim')(trim)
    ConfigMap(_func='collapse')(collapse)
    for argmap, params in loadConfig(sys.argv, (ProDuSe, trim, collapse), batchExpression='samples'):
        global config
        config = argmap
        ProDuSe(**argmap[ProDuSe])
        break #TODO testing only, delete