#!/usr/bin/env python3
import sys, os, tempfile, shutil, gzip

from collapse import collapse
from clip import clip
from trim import trim

from multiprocessing import Process
from subprocess import Popen, PIPE

from configutator import ConfigMap, loadConfig
from valitator import Path, UnsignedInt, PathOrNone, Executable

@ConfigMap(fastq1='fastqs[0]', fastq2='fastqs[1]', bwaArgs='bwa', output='join(``, [name, `.bam`])')
def ProDuSe(fastq1: Path, fastq2: Path, reference: str, bed: PathOrNone, output: str, bwaArgs: Executable):
    """

    :param fastq1: Path to the first fastq
    :param fastq2: Path to the second fastq
    :param reference: Path to the reference fasta
    :param bed: Path to a bed file containing regions to restrict reads to (optional)
    :param output: Path to output file
    :param bwaArgs:
    :return:
    """
    global config

    temp_dir = tempfile.mkdtemp()
    sys.stderr.write("Pipes created in {}\n".format(temp_dir))
    trimOut1Path = os.path.join(temp_dir, 'trimOut1')
    os.mkfifo(trimOut1Path)
    trimOut2Path = os.path.join(temp_dir, 'trimOut2')
    os.mkfifo(trimOut2Path)

    # Check if fastqs compressed
    inFile1 = open(fastq1, 'rb')
    if inFile1.peek(2)[:2] == b'\037\213':
        inFile1 = gzip.open(inFile1)
    inFile2 = open(fastq2, 'rb')
    if inFile2.peek(2)[:2] == b'\037\213':
        inFile2 = gzip.open(inFile2)

    # Start trim subprocesses
    trimProc1 = Process(target=trim, args=(inFile1, open(trimOut1Path, 'wb')))
    trimProc1.start()
    trimProc2 = Process(target=trim, args=(inFile2, open(trimOut2Path, 'wb')))
    trimProc2.start()

    # Start bwa subprocess
    if isinstance(bwaArgs, str):
        bwaArgs = bwaArgs.split(' ')
    sys.stderr.write("Starting BWA subprocess..\n")
    bwa = Popen(bwaArgs + ['-C', reference, trimOut1Path, trimOut2Path], stdout=PIPE)
    clipIn = bwa.stdout

    if bed:
        #Start samtools view to restrict to coordinates
        sys.stderr.write("Starting samtools view subprocess..\n")
        view = Popen(["samtools", "view", "-uL", bed, "-"], stdin=clipIn, stdout=PIPE)
        clipIn = view.stdout

    # Start sort by coord
    sys.stderr.write("Starting sort subprocess..\n")
    sort = Popen(["samtools", "sort", "-l0", "-"], stdout=PIPE, stdin=PIPE)

    #Start clipping subprocess
    sys.stderr.write("Starting clipping subprocess..\n")
    clipProc = Process(target=clip, args=(clipIn, sort.stdin))
    clipProc.start()

    sys.stderr.write("Starting collapse subprocess..\n")
    collapseProc = Process(target=collapse, kwargs=dict(inStream=sort.stdout, outStream=open(output, 'wb+'), **config[collapse]))
    collapseProc.start()

    clipProc.join()
    collapseProc.join()
    bwa.wait()
    inFile1.close()
    inFile2.close()
    shutil.rmtree(temp_dir)

if __name__ == "__main__":
    ConfigMap(_func='trim')(trim)
    ConfigMap(_func='collapse')(collapse)
    for argmap in loadConfig(sys.argv, (ProDuSe, trim, collapse)):
        global config
        config = argmap
        ProDuSe(**argmap[ProDuSe])
        break #TODO testing only, delete