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
  run_produse     Runs all steps of the ProDuSe pipeline listed below
  trim            Removes barcodes from FASTQ files
  clip            Clips bases which overlap in both the forward and reverse reads in paired-end sequencing
  collapse        Identifies duplicates,and merges duplicates into a consensus
  call            Identifies variants
  All commands accept -h for their specific usage.

  
The relevant publications can be found here:

- 
- 

Collapse
--------

Collapse reduces records with the same reference start position, forward/reverse flag and barcode within a specified Hamming distance to a single concensus family record. The families are first collapsed if the barcodes are identical into "proto-families". These proto-families barcodes are then compared to each other producing a graph with node weights of the proto-family size and edge weights of their Hamming distances. The graph is reduced to a minimum spanning forest where any edges above threshold are removed. Each tree in the forest is further reduced to a diameter twice the threshold while maximising each subtrees node weight. 
A concensus of each tree is then finalised.

Each family record is output with two additional tags:

- fQ: Integer array containing Phred score of wrong base chosen during collapse
- fC: Integer array containing pairs (simmilar to a CIGAR) of count and depth representing, from the start of the alignment, the depth of the family at a position

Call
----

