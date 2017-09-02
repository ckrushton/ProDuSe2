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
  adapter_predict Identifies the degenerate barcode used in a set of FASTQ files
  run             Runs all steps of the ProDuSe pipeline listed below
  trim            Removes barcodes from FASTQ files
  collapse        Identifies duplicates,and merges duplicates into a consensus
  call            Identifies variants
  All commands accept -h for their specific usage.

  
The relevant publications can be found here:

- 
- 

Note:
-----

- ProDuSe makes use of `ClipOverlap <https://github.com/innovate-invent/clip>`_ in its pipeline with arguments -abt8.
- bwa mem is the aligner used.

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

Configuration
-------------

While all arguments can be passed via command line, the main pipeline requires so many that a config file is recommended. Config files are in YAML format. See `configutator <https://github.com/innovate-invent/configutator>`_ for more details.

Here is an example config.yaml with all parameters specified::

  default: &default
    output: ./
    reference: Ref.fa
    bwa: bwa mem -t 10
    samtools: samtools

    trim: &trim_default
      barcode_distance: 3

    collapse: &collapse_default
      barcode_mask: "SSS1111111111110"
      barcode_distance: 3

    filter: &filter_default
      min_molecules: 2
      mutant_molecules: 3
      variant_allele_fraction_threshold: 0.01
      min_reads_per_uid: 1

  samples:
    - name: test1
      <<: *default     # Inherit default parameters
      fastqs:
        - test/test1/onm_R1.fastq
        - test/test1/onm_R1.fastq
      trim:
        <<: *trim_default
        barcode_sequence: NNNWSMRWSYWKMWWT

    - name: test2
      <<: *default     # Inherit default parameters
      fastqs:
        - test/test2/afnm_R1.fastq
        - test/test2/afnm_R2.fastq
      trim:
        <<: *trim_default
        barcode_sequence: NNNWSMRWSYWKMWWT
