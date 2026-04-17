#!/usr/bin/env python3
"""
metablastbin.py — Multi-sample taxonomic binning from BLAST tabular output

Processes a folder of BLAST TSV files. For each contig in each sample it
picks the best BLAST hit (highest bitscore) and records the assigned
organism. All results are written to a single combined summary TSV.

Optionally, per-sample FASTA files can also be provided in the same folder
(matched by basename) so that unassigned contigs are captured and per-bin
FASTA files can be written.

Folder layout expected
----------------------
    blast_dir/
        sample1.tsv          ← BLAST output
        sample1.fasta        ← assembly (optional, matched by basename)
        sample2.tsv
        sample2.fasta
        ...

FASTA extensions recognised: .fa  .fasta  .fna

Column layout of BLAST TSV (0-based)
-------------------------------------
  0   qseqid    contig name
  1   qlen      contig length
  2   sseqid    subject accession
  3   slen      subject length
  4   pident    % identity
  5   length    alignment length
  6   evalue
  7   qstart / 8 qend / 9 sstart / 10 send / 11 N/A
  12  stitle    full subject title  ← organism name derived from here
  13  staxid
  14  bitscore

Usage
-----
    python metablastbin.py -d blast_dir/ -o results/ [options]

Options
-------
  -d / --blast-dir   Folder containing BLAST TSV files (required)
  -o / --outdir      Output folder (default: binner_out/)
  -e / --evalue      Max e-value  (default: 1e-5)
  -p / --pident      Min % identity (default: 80.0)
  -a / --alen        Min alignment length in bp (default: 0)
  --tax-level        species | genus  (default: species)
  --write-fasta      Write per-bin FASTA files (needs matching assembly FASTAs)
  --print-summary    Print the combined summary table to stdout as well
  --sep              Column separator in BLAST files (default: TAB)

Output
------
  results/
      combined_summary.tsv          ← one row per contig × sample
      <sample>/                     ← sub-folder per sample (only if --write-fasta)
          <organism>.fa
          unassigned.fa
"""

import argparse
import os
import re
import sys
from collections import defaultdict

# ---------------------------------------------------------------------------
# BLAST column indices
# ---------------------------------------------------------------------------
COL_QSEQID   = 0
COL_QLEN     = 1
COL_SSEQID   = 2
COL_PIDENT   = 4
COL_ALEN     = 5
COL_EVALUE   = 6
COL_STITLE   = 12
COL_BITSCORE = 14

FASTA_EXTS = (".fa", ".fasta", ".fna")
BLAST_EXTS = (".tsv", ".txt", ".blast", ".blastout", ".out")


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Multi-sample BLAST-based contig binner.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("-d", "--blast-dir",   required=True,
                   help="Folder containing BLAST TSV files")
    p.add_argument("-o", "--outdir",      default="binner_out",
                   help="Output folder (default: binner_out/)")
    p.add_argument("-e", "--evalue",      type=float, default=1e-5,
                   help="Max e-value (default: 1e-5)")
    p.add_argument("-p", "--pident",      type=float, default=80.0,
                   help="Min %% identity (default: 80.0)")
    p.add_argument("-a", "--alen",        type=int,   default=0,
                   help="Min alignment length in bp (default: 0)")
    p.add_argument("--tax-level",         default="species",
                   choices=["species", "genus"],
                   help="Taxonomic level for binning (default: species)")
    p.add_argument("--fasta-dir",          default=None,
                   help="Folder containing assembly FASTAs if different from --blast-dir. "
                        "FASTAs are matched to BLAST files by basename "
                        "(e.g. sample1.tsv ↔ sample1.fasta). "
                        "If omitted, FASTAs are looked up in --blast-dir.")
    p.add_argument("--write-fasta",       action="store_true",
                   help="Write per-bin FASTA files. "
                        "Requires a matching assembly FASTA for every BLAST file.")
    p.add_argument("--print-summary",     action="store_true",
                   help="Print combined summary table to stdout")
    p.add_argument("--sep",               default="\t",
                   help="Column separator in BLAST files (default: TAB)")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Organism name helpers
# ---------------------------------------------------------------------------

def clean_title(title: str) -> str:
    title = re.sub(r"^(MAG\s+)?TPA[_\w]*:\s*", "", title, flags=re.I)
    title = re.sub(r"^UNVERIFIED[_\w]*:\s*", "", title, flags=re.I)
    return title.strip()


def extract_organism(title: str, level: str) -> str:
    title = clean_title(title)
    stop_re = re.compile(
        r",|\bstrain\b|\bisolate\b|\bchromosome\b|\bplasmid\b|\bscaffold\b|\bcontig\b",
        re.I,
    )
    m = stop_re.search(title)
    trimmed = title[: m.start()].strip() if m else title
    words = trimmed.split()

    if level == "genus":
        name = words[0] if words else "Unknown"
    else:  # species
        name = " ".join(words[:2]) if len(words) >= 2 else (words[0] if words else "Unknown")

    name = re.sub(r"[^\w\s.\-]", "", name).strip()
    name = re.sub(r"\s+", "_", name)
    return name or "Unknown"


def safe_filename(name: str) -> str:
    return re.sub(r"[^\w.\-]", "_", name)


# ---------------------------------------------------------------------------
# FASTA parsing
# ---------------------------------------------------------------------------

def parse_fasta(path: str) -> dict:
    """Return { contig_id: {"seq": str, "length": int} }"""
    contigs = {}
    current_id = None
    chunks = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_id is not None:
                    seq = "".join(chunks)
                    contigs[current_id] = {"seq": seq, "length": len(seq)}
                current_id = line[1:].split()[0]
                chunks = []
            elif line:
                chunks.append(line)
    if current_id is not None:
        seq = "".join(chunks)
        contigs[current_id] = {"seq": seq, "length": len(seq)}
    return contigs


def write_fasta_bin(path: str, entries: list):
    """Write a FASTA file. entries = [(header, seq), ...]"""
    with open(path, "w") as fh:
        for header, seq in entries:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")


# ---------------------------------------------------------------------------
# BLAST parsing
# ---------------------------------------------------------------------------

def parse_blast(path: str, sep: str, max_evalue: float,
                min_pident: float, min_alen: int) -> tuple:
    """
    Returns:
        best_hit  { contig_id: { qlen, evalue, pident, alen, sseqid, stitle, bitscore } }
        stats     { n_lines, n_passed, n_skipped }
    """
    best_hit = {}
    n_lines = n_passed = n_skipped = 0

    required_cols = max(COL_QSEQID, COL_QLEN, COL_SSEQID,
                        COL_PIDENT, COL_ALEN, COL_EVALUE, COL_STITLE) + 1

    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            n_lines += 1
            cols = line.split(sep)

            if len(cols) < required_cols:
                n_skipped += 1
                continue

            try:
                qseqid   = cols[COL_QSEQID].strip()
                qlen     = int(cols[COL_QLEN].strip())
                sseqid   = cols[COL_SSEQID].strip()
                pident   = float(cols[COL_PIDENT].strip())
                alen     = int(cols[COL_ALEN].strip())
                evalue   = float(cols[COL_EVALUE].strip())
                stitle   = cols[COL_STITLE].strip()
                bitscore = float(cols[COL_BITSCORE].strip()) if len(cols) > COL_BITSCORE else 0.0
            except (ValueError, IndexError):
                n_skipped += 1
                continue

            # Filter
            if evalue > max_evalue or pident < min_pident or alen < min_alen:
                continue

            n_passed += 1

            # Keep best hit per contig
            prev = best_hit.get(qseqid)
            if prev is None or bitscore > prev["bitscore"] or (
                bitscore == prev["bitscore"] and evalue < prev["evalue"]
            ):
                best_hit[qseqid] = dict(
                    qlen=qlen, evalue=evalue, pident=pident,
                    alen=alen, sseqid=sseqid, stitle=stitle, bitscore=bitscore,
                )

    return best_hit, {"n_lines": n_lines, "n_passed": n_passed, "n_skipped": n_skipped}


# ---------------------------------------------------------------------------
# File discovery
# ---------------------------------------------------------------------------

def find_blast_files(blast_dir: str) -> list:
    """Return sorted list of BLAST TSV paths found in blast_dir."""
    files = []
    for fname in sorted(os.listdir(blast_dir)):
        ext = os.path.splitext(fname)[1].lower()
        if ext in BLAST_EXTS:
            files.append(os.path.join(blast_dir, fname))
    return files


def find_fasta_for(blast_path: str, fasta_dir: str = None) -> str | None:
    """
    Look for a FASTA with the same basename as the BLAST file.
    Search order:
      1. fasta_dir  (if provided)
      2. same directory as the BLAST file
    """
    basename = os.path.splitext(os.path.basename(blast_path))[0]
    search_dirs = []
    if fasta_dir:
        search_dirs.append(fasta_dir)
    search_dirs.append(os.path.dirname(blast_path))

    for d in search_dirs:
        for ext in FASTA_EXTS:
            candidate = os.path.join(d, basename + ext)
            if os.path.isfile(candidate):
                return candidate
    return None


def sample_name_from_path(blast_path: str) -> str:
    """e.g. /data/sample1.tsv  →  sample1"""
    return os.path.splitext(os.path.basename(blast_path))[0]


# ---------------------------------------------------------------------------
# Per-sample processing
# ---------------------------------------------------------------------------

def process_sample(blast_path: str, args) -> list:
    """
    Process one sample. Returns a list of row dicts for the combined summary.
    Also writes per-bin FASTAs if --write-fasta is set.
    """
    sample = sample_name_from_path(blast_path)
    print(f"[MetaBlastBin]  Sample: {sample}", file=sys.stderr)

    # --- BLAST ---
    best_hit, stats = parse_blast(
        blast_path, args.sep, args.evalue, args.pident, args.alen
    )
    print(
        f"[MetaBlastBin]    lines={stats['n_lines']:,}  "
        f"passed={stats['n_passed']:,}  "
        f"contigs_with_hit={len(best_hit):,}  "
        f"skipped={stats['n_skipped']:,}",
        file=sys.stderr,
    )

    # --- Optional FASTA ---
    fasta_path = find_fasta_for(blast_path, fasta_dir=getattr(args, "fasta_dir", None))
    fasta_contigs = {}
    if fasta_path:
        fasta_contigs = parse_fasta(fasta_path)
        print(f"[MetaBlastBin]    FASTA: {os.path.basename(fasta_path)} "
              f"({len(fasta_contigs):,} contigs)", file=sys.stderr)
    elif args.write_fasta:
        # FASTA is required for MAG extraction — hard stop
        print(
            f"\n[MetaBlastBin] ERROR: --write-fasta is set but no matching FASTA found for:\n"
            f"  BLAST file : {blast_path}\n"
            f"  Looked for : {os.path.splitext(os.path.basename(blast_path))[0]}"
            f"{{.fa,.fasta,.fna}} in:\n"
            f"    {getattr(args, 'fasta_dir', None) or '(not set)'}\n"
            f"    {os.path.dirname(blast_path)}\n"
            f"\nFix: make sure your FASTA basename matches the BLAST basename, e.g.:\n"
            f"  sample1.tsv  ↔  sample1.fasta\n"
            f"Or use --fasta-dir to point to the folder containing your FASTAs.",
            file=sys.stderr,
        )
        sys.exit(1)

    # --- Build summary rows ---
    rows = []

    for contig_id, hit in best_hit.items():
        organism = extract_organism(hit["stitle"], args.tax_level)
        rows.append(dict(
            sample      = sample,
            contig_id   = contig_id,
            contig_len  = hit["qlen"],
            organism    = organism,
            full_title  = hit["stitle"],
            accession   = hit["sseqid"],
            pident      = round(hit["pident"], 2),
            evalue      = hit["evalue"],
            bitscore    = hit["bitscore"],
            status      = "assigned",
        ))

    # Unassigned (only if FASTA provided)
    unassigned_ids = []
    if fasta_contigs:
        unassigned_ids = [cid for cid in fasta_contigs if cid not in best_hit]
        for cid in unassigned_ids:
            rows.append(dict(
                sample      = sample,
                contig_id   = cid,
                contig_len  = fasta_contigs[cid]["length"],
                organism    = "unassigned",
                full_title  = "",
                accession   = "",
                pident      = "",
                evalue      = "",
                bitscore    = "",
                status      = "unassigned",
            ))
        print(f"[MetaBlastBin]    unassigned={len(unassigned_ids):,}", file=sys.stderr)

    # --- Optional per-bin FASTAs ---
    if args.write_fasta and fasta_contigs:
        sample_outdir = os.path.join(args.outdir, safe_filename(sample))
        os.makedirs(sample_outdir, exist_ok=True)

        # Group assigned contigs by organism
        bins = defaultdict(list)
        for row in rows:
            if row["status"] == "assigned":
                bins[row["organism"]].append(row["contig_id"])

        for org, contig_ids in bins.items():
            entries = []
            for cid in sorted(contig_ids):
                entry = fasta_contigs.get(cid)
                if entry:
                    entries.append((f"{cid} {org}", entry["seq"]))
            if entries:
                fa_path = os.path.join(sample_outdir, safe_filename(org) + ".fa")
                write_fasta_bin(fa_path, entries)

        # Unassigned FASTA
        if unassigned_ids:
            ua_entries = []
            for cid in sorted(unassigned_ids):
                entry = fasta_contigs.get(cid)
                if entry:
                    ua_entries.append((cid, entry["seq"]))
            if ua_entries:
                write_fasta_bin(
                    os.path.join(sample_outdir, "unassigned.fa"), ua_entries
                )

        print(f"[MetaBlastBin]    FASTA bins → {sample_outdir}/", file=sys.stderr)

    return rows


# ---------------------------------------------------------------------------
# Combined summary output
# ---------------------------------------------------------------------------

SUMMARY_COLS = [
    "sample", "contig_id", "contig_len", "organism",
    "accession", "pident", "evalue", "bitscore", "full_title", "status",
]


def write_combined_summary(all_rows: list, outdir: str) -> str:
    os.makedirs(outdir, exist_ok=True)
    out_path = os.path.join(outdir, "combined_summary.tsv")
    with open(out_path, "w") as fh:
        fh.write("\t".join(SUMMARY_COLS) + "\n")
        for row in all_rows:
            fh.write("\t".join(str(row.get(c, "")) for c in SUMMARY_COLS) + "\n")
    return out_path


def print_summary_table(all_rows: list):
    """Print a compact per-sample × organism count table to stdout."""
    # Aggregate: sample → organism → contig count
    agg = defaultdict(lambda: defaultdict(int))
    for row in all_rows:
        agg[row["sample"]][row["organism"]] += 1

    # Collect all organisms across all samples
    all_orgs = sorted({row["organism"] for row in all_rows})
    samples   = sorted(agg.keys())

    org_w  = max(len(o) for o in all_orgs) + 2
    org_w  = max(org_w, 20)
    smp_w  = max(len(s) for s in samples) + 2
    smp_w  = max(smp_w, 12)
    col_w  = 8

    # Header
    header = f"{'Organism':<{org_w}}" + "".join(f"{s:>{smp_w}}" for s in samples)
    print("\n" + "=" * len(header))
    print(header)
    print("-" * len(header))
    for org in all_orgs:
        row_str = f"{org:<{org_w}}"
        for s in samples:
            cnt = agg[s].get(org, 0)
            row_str += f"{cnt:>{smp_w},}"
        print(row_str)
    print("=" * len(header))
    total = len(all_rows)
    print(f"Total rows: {total:,}\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    blast_files = find_blast_files(args.blast_dir)
    if not blast_files:
        print(f"[MetaBlastBin] ERROR: No BLAST files found in {args.blast_dir}", file=sys.stderr)
        print(f"[MetaBlastBin] Recognised extensions: {BLAST_EXTS}", file=sys.stderr)
        sys.exit(1)

    print(f"[MetaBlastBin] Found {len(blast_files)} BLAST file(s) in {args.blast_dir}",
          file=sys.stderr)
    print(f"[MetaBlastBin] Filters  — e-value ≤ {args.evalue}  "
          f"identity ≥ {args.pident}%  alignment ≥ {args.alen} bp", file=sys.stderr)
    print(f"[MetaBlastBin] Tax level: {args.tax_level}", file=sys.stderr)
    print(f"[MetaBlastBin] Output  → {args.outdir}/", file=sys.stderr)
    print(file=sys.stderr)

    all_rows = []
    for blast_path in blast_files:
        rows = process_sample(blast_path, args)
        all_rows.extend(rows)
        print(file=sys.stderr)

    summary_path = write_combined_summary(all_rows, args.outdir)
    print(f"[MetaBlastBin] Combined summary → {summary_path}", file=sys.stderr)
    print(f"[MetaBlastBin] Total rows written: {len(all_rows):,}", file=sys.stderr)

    if args.print_summary:
        print_summary_table(all_rows)


if __name__ == "__main__":
    main()
