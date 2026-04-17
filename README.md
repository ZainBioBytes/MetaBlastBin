# MetaBlastBin

**BLAST-based metagenomic contig binner** — assigns contigs to organisms from BLAST tabular output and extracts per-organism MAGs across multiple samples.

---

## What it does

MetaBlastBin takes your BLAST results and assembly FASTAs and bins every contig into the organism it most likely came from. For each contig it picks the single best BLAST hit (highest bitscore), extracts the organism name from the subject title, and groups all contigs assigned to the same organism into one FASTA file — effectively producing metagenome-assembled genomes (MAGs) directly from BLAST identity rather than coverage-based methods.

It processes an entire folder of samples in one run and produces a single combined summary TSV across all samples.

```
blast_dir/
    sample1.tsv  +  sample1.fasta  →  bins/sample1/Comamonas_kerstersii.fa
                                       bins/sample1/Staphylococcus_chromogenes.fa
                                       bins/sample1/unassigned.fa
    sample2.tsv  +  sample2.fasta  →  bins/sample2/Escherichia_coli.fa
                                       ...
                                   →  bins/combined_summary.tsv
```

---

## Requirements

- Python 3.9+
- No external dependencies — standard library only

---

## Installation

```bash
git clone https://github.com/ZainBioBytes/MetaBlastBin.git
cd MetaBlastBin
python metablastbin.py --help
```

---

## BLAST format expected

MetaBlastBin expects BLAST tabular output with the following custom columns (produced with `-outfmt`):

```
-outfmt "6 qseqid qlen sseqid slen pident length evalue qstart qend sstart send scomnames stitle staxid bitscore"
```

| Col | Field | Description |
|-----|-------|-------------|
| 0 | qseqid | Contig name |
| 1 | qlen | Contig length |
| 2 | sseqid | Subject accession |
| 3 | slen | Subject length |
| 4 | pident | % identity |
| 5 | length | Alignment length |
| 6 | evalue | E-value |
| 7–11 | qstart qend sstart send (extra) | Coordinates |
| 12 | stitle | Full subject title ← organism name extracted from here |
| 13 | staxid | NCBI taxon ID |
| 14 | bitscore | Bitscore |

---

## Usage

### Basic — summary only (no FASTA needed)

```bash
python metablastbin.py \
    -d blast_results/ \
    -o output/
```

### Full — extract MAG FASTAs

FASTAs and BLAST files in the **same folder** (matched by basename):

```bash
python metablastbin.py \
    -d blast_results/ \
    -o output/ \
    --write-fasta
```

FASTAs in a **separate folder**:

```bash
python metablastbin.py \
    -d blast_results/ \
    --fasta-dir assemblies/ \
    -o output/ \
    --write-fasta
```

Print a cross-sample organism count table to the terminal:

```bash
python metablastbin.py \
    -d blast_results/ \
    --fasta-dir assemblies/ \
    -o output/ \
    --write-fasta \
    --print-summary
```

---

## All options

| Flag | Default | Description |
|------|---------|-------------|
| `-d` / `--blast-dir` | required | Folder containing BLAST TSV files |
| `-o` / `--outdir` | `binner_out/` | Output folder |
| `--fasta-dir` | same as `--blast-dir` | Folder with assembly FASTAs (if separate) |
| `-e` / `--evalue` | `1e-5` | Maximum e-value to accept a hit |
| `-p` / `--pident` | `80.0` | Minimum % identity |
| `-a` / `--alen` | `0` | Minimum alignment length (bp) |
| `--tax-level` | `species` | Bin at `species` (2 words) or `genus` (1 word) level |
| `--write-fasta` | off | Write per-bin FASTA files (MAGs) |
| `--print-summary` | off | Print organism × sample count table to stdout |
| `--sep` | TAB | Column separator in BLAST files |

---

## File matching (BLAST ↔ FASTA)

The script matches each BLAST file to its assembly FASTA by **basename**:

```
sample1.tsv   ↔   sample1.fa  /  sample1.fasta  /  sample1.fna
calf_D01.blast ↔  calf_D01.fasta
```

The name before the extension must be **identical** (case-sensitive). If a FASTA is not found when `--write-fasta` is set, the script will exit with a clear error message.

---

## Output

### `combined_summary.tsv`

One row per contig per sample:

| sample | contig_id | contig_len | organism | accession | pident | evalue | bitscore | full_title | status |
|--------|-----------|------------|----------|-----------|--------|--------|----------|------------|--------|
| sample1 | contig_1 | 86801 | Comamonas_kerstersii | CP020121 | 98.18 | 0.0 | 28.0 | Comamonas kerstersii strain... | assigned |
| sample1 | contig_orphan | 11 | unassigned | | | | | | unassigned |

### Per-bin FASTAs (with `--write-fasta`)

```
output/
    combined_summary.tsv
    sample1/
        Comamonas_kerstersii.fa
        Staphylococcus_chromogenes.fa
        Caudoviricetes_sp..fa
        unassigned.fa          ← contigs with no passing BLAST hit
    sample2/
        Escherichia_coli.fa
        Klebsiella_pneumoniae.fa
        Comamonas_kerstersii.fa
        unassigned.fa
```

---

## Try the example

```bash
python metablastbin.py \
    -d usage_example/blast_input/ \
    --fasta-dir usage_example/assemblies/ \
    -o usage_example/my_output/ \
    --write-fasta \
    --print-summary
```

Expected output is in `usage_example/expected_output/`.

---

## How binning works

1. All HSPs for a contig are read from the BLAST file
2. The single best hit is selected (highest bitscore; tie-break: lowest e-value)
3. Hits that fail the e-value, identity, or alignment length filters are discarded
4. The organism name is extracted from the `stitle` column — prefixes like `MAG TPA_asm:` are stripped, and suffixes like `strain 8943, complete genome` are dropped, leaving a clean species name
5. All contigs assigned to the same organism are grouped into one bin
6. If an assembly FASTA is provided, contigs with no passing hit go to `unassigned.fa`

> **Note:** This is winner-takes-all binning. A contig is assigned entirely to one organism based on its top hit. Chimeric contigs are not split.

---

## Citation / acknowledgements

If you use MetaBlastBin in your work, please cite this repository.

---

## License

MIT License — see [LICENSE](LICENSE).
