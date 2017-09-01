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

