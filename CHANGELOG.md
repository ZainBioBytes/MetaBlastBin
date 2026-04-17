# Changelog

All notable changes to MetaBlastBin will be documented here.

## [1.0.0] - 2025-04-17

### Initial release

- Multi-sample batch processing from a folder of BLAST TSV files
- Best-hit selection per contig (highest bitscore, tie-break by e-value)
- Organism name extraction from `stitle` with prefix/suffix stripping
- Species and genus level binning via `--tax-level`
- Per-bin MAG FASTA extraction with `--write-fasta`
- Support for FASTAs in a separate folder via `--fasta-dir`
- Single combined summary TSV across all samples
- Unassigned contig tracking when assembly FASTA is provided
- Clear error messages when FASTA files are missing
- No external dependencies — Python standard library only
