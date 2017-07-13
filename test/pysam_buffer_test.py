import pysam, multiprocessing, sys

def p(pR):
    c = pysam.AlignmentFile(pR.fileno())
    d = pysam.AlignmentFile(sys.stdout, 'w', template=c)
    for r in c.fetch(until_eof=True): d.write(r)
    c.close()
    d.close()

a = pysam.AlignmentFile('/home/ncm3/data/singlePair250.bam')
pR, pW = multiprocessing.Pipe(False)
b = pysam.AlignmentFile(pW.fileno(), 'w', template=a)

proc = multiprocessing.Process(target=p, args=(pR,))
proc.start()

for r in a: b.write(r)
a.close()
b.close()
pW.close()

#proc = multiprocessing.Process(target=p, args=(pR,))
#proc.start()

#p(pR)

pass