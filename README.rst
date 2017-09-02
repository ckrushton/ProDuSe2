ProDuSe2
--------

Variant caller that utilises semi-degenerate barcoded adapter libraries.

To install: ``pip install produse``
Or pull from the official repo::

  git clone git@github.com:morinlab/ProDuSe.git
  cd ProDuSe
  python3 setup.py install
  
To use:
-------

::

  $ produse
  Usage:   produse <command> [options]
  Commands:
  adapter_predict Identifies the degerate barcode used in a set of FASTQ files
  run             Runs all steps of the ProDuSe pipeline listed below
  trim            Removes barcodes from FASTQ files
  collapse        Identifies duplicates,and merges duplicates into a consensus
  call            Identifies variants
  All commands accept -h for their specific usage.

  
The relevant publications can be found here:

- 
- 

Note: ProDuSe makes use of `ClipOverlap <https://github.com/innovate-invent/clip>`_ in its pipeline with parameters -abt8

Collapse
--------

Collapse reduces records with the same reference start position, forward/reverse flag and barcode within a specified Hamming distance to a single concensus family record. A consensus family record is essentially a pileup of the family member records and the most common operation/base at each position is retained.
The families are first collapsed if the barcodes are identical into "proto-families". These proto-families barcodes are then compared to each other producing a graph with node weights of the proto-family size and edge weights of their Hamming distances. The graph is reduced to a minimum spanning forest where any edges above threshold are removed. Each tree in the forest is further reduced to a diameter twice the threshold while maximising each subtrees node weight. 
A concensus of each tree is then finalised.

Each family record is output with two additional tags:

- fQ: Integer array containing Phred score of wrong base chosen during collapse
- fC: Integer array containing pairs (simmilar to a CIGAR) of count and depth representing, from the start of the alignment, the depth of the family at a position

Call
----

