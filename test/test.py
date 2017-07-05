def test_condense():
    from importlib import reload
    import pysam
    import condense

    inf = pysam.AlignmentFile("/home/ncm3/data/normal.bam", "rb")
    outf = pysam.AlignmentFile("/home/ncm3/data/normal_condensed.bam", "wb", inf)

    condense.condense(inf, outf)

    inf.close()
    outf.close()

def test_configutator():
    import configutator
    import produse

    configutator.generate([produse.ProDuSe])

test_configutator()