# gfflu changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.0.2](https://github.com/CFIA-NCFAD/gfflu/releases/tag/0.0.2) - 2023-08-09

### Fixed

* Handling of CDS with joins (e.g. `join(24..49,738..1005)` for M2 of MH085255) is now properly combined and translated to amino acid sequence (#3)

### Added

* Amino acid FASTA output for all CDS and other features that can be translated to peptides (`{outdir}/{prefix}.faa`)
* `--prefix/-p` CLI option for specifying output file prefix

### Changed

* `peptides.py`: convert ha_signal_peptides and pb1_f2_peptides to tuples instead of list; should not be modified during runtime

## [0.0.1](https://github.com/CFIA-NCFAD/gfflu/releases/tag/0.0.1) - 2023-06-09

First release of gfflu.
