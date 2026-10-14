## v2.1.0 - 2026-09-21

## Updates

- Alter output CS file setting CS tag.
- Allow separation of sequences generated.
- File with alleles BLAST

## Fixed

- Correction of bornes temp
- Correction of alleles blast


## v2.0.2 - 2026-09-15

### Fixed

- Clippy correction
- Fix bornes and locus delimitation
- Both haplotype of a same locus cannot be on the same chromosome anymore

## v2.0.1 - 2026-08-20

### Fixed

- Plus/Minus strand for gene
- Update of MD parsing (update of extended-htslib)

### Soon (version 2.1 )

- Alter output CS file setting CS tag.

## v2.0.0 - 2026-08-19

IMGT/StatAssembly v2.0.0 is a major version update bringing new functionalities.

### Added
- Add automatic detection of locus and bornes
- Add automatic generation of a gene list
- Add detection of soft clips and mean coverage.
- Automatic acceptation or rejection of locus and genes based on IMGT criterias.
- Creation of a small BAM file containing BAM alignment for all loci, as well as a compressed fasta file of all the sequences and the list of validated alleles.

### Changed
- IMGT criterias to accept or reject a locus and/or an allele.
