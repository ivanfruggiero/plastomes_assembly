#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path


def read_fasta(path: Path):
    """Read a FASTA file and return a list of (name, sequence) tuples."""
    recs = []
    name = None
    seq = []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(seq).upper()))
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if name is not None:
            recs.append((name, "".join(seq).upper()))
    return recs


def write_fasta(path: Path, header: str, seq: str, width: int = 80):
    """Write a single-record FASTA file with fixed line width."""
    with path.open("w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(seq), width):
            f.write(seq[i : i + width] + "\n")


def revcomp(seq: str) -> str:
    """Reverse-complement a DNA sequence (supports A,C,G,T,N)."""
    tr = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(tr)[::-1]


def build_ref_full(ref_records) -> str:
    """
    Build a full reference sequence from FASTA records.

    Accepted formats:
      1) single-record FASTA (most common)
      2) quadripartite FASTA with 4 records (e.g., LSC / IRb / SSC / IRa)

    In the quadripartite case, the full reference is obtained by concatenating
    records in the order they appear in the FASTA file. This implies that the
    reference should be exported/created in the desired LSC-IRb-SSC-IRa order.
    """
    if not ref_records:
        raise ValueError("Reference FASTA has no records.")
    return "".join(seq for _, seq in ref_records)


def rotate(seq: str, idx: int) -> str:
    """Rotate sequence so that seq[idx] becomes the first base."""
    return seq[idx:] + seq[:idx]


def try_anchor(seq: str, ref_full: str, anchor_len: int):
    """
    Try to find the reference start anchor within seq (forward or reverse-complement).
    Returns (normalized_seq, orientation, used_anchor_len, anchor_pos) or (None, ...).
    """
    anchor = ref_full[:anchor_len]
    i = seq.find(anchor)
    if i != -1:
        return rotate(seq, i), "forward", anchor_len, i

    rc = revcomp(seq)
    j = rc.find(anchor)
    if j != -1:
        return rotate(rc, j), "revcomp", anchor_len, j

    return None, "not_found", anchor_len, -1


def normalize_to_reference(seq: str, ref_full: str, anchor_len: int, min_anchor: int, step: int):
    """
    Normalize a plastome by rotating it to match the reference start (anchor-based).

    Strategy:
      - Take an anchor from the start of the reference (length = anchor_len)
      - Search it in the sample sequence; if not found, try reverse-complement
      - If still not found, decrease anchor length by 'step' until min_anchor

    Returns (normalized_seq, orientation, used_anchor_len, anchor_pos) or (None, ...)
    """
    L = anchor_len
    while L >= min_anchor:
        out, orient, usedL, pos = try_anchor(seq, ref_full, L)
        if out is not None:
            return out, orient, usedL, pos
        L -= step
    return None, "not_found", 0, -1


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Normalize plastome FASTA sequences by rotating them to match the reference start (anchor-based). "
            "Reference can be either a single-record FASTA or a quadripartite FASTA with 4 records "
            "(e.g., LSC/IRb/SSC/IRa), which will be concatenated in file order."
        )
    )
    ap.add_argument("--ref", required=True, help="Reference FASTA (single record or 4-record quadripartite).")
    ap.add_argument("--in_dir", required=True, help="Input directory containing sample FASTA files.")
    ap.add_argument("--out_dir", required=True, help="Output directory for normalized FASTA files.")
    ap.add_argument("--pattern", default="*.fasta", help="Glob pattern for input FASTA (default: *.fasta).")
    ap.add_argument("--anchor_len", type=int, default=200, help="Anchor length from reference start (default: 200).")
    ap.add_argument("--min_anchor", type=int, default=75, help="Minimum anchor length (default: 75).")
    ap.add_argument("--step", type=int, default=25, help="Decrease anchor by this step (default: 25).")
    args = ap.parse_args()

    ref_path = Path(args.ref)
    in_dir = Path(args.in_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if not ref_path.exists():
        raise SystemExit(f"Reference FASTA not found: {ref_path}")
    if not in_dir.exists():
        raise SystemExit(f"Input directory not found: {in_dir}")

    ref_records = read_fasta(ref_path)
    ref_full = build_ref_full(ref_records)

    if args.anchor_len > len(ref_full):
        raise SystemExit(f"--anchor_len ({args.anchor_len}) is longer than reference length ({len(ref_full)}).")

    in_files = sorted(in_dir.glob(args.pattern))
    if not in_files:
        raise SystemExit(f"No input FASTA files found in {in_dir} with pattern {args.pattern}")

    report = ["file\tsample\tstatus\torientation\tanchor_len\tanchor_pos\tseq_len_bp\tout_file"]

    for fp in in_files:
        sample = fp.stem
        recs = read_fasta(fp)
        if not recs:
            report.append(f"{fp.name}\t{sample}\tempty\t-\t0\t-1\t0\t-")
            continue

        # Use the first record (typical GetOrganelle FASTA output is single-record)
        seq = recs[0][1]
        norm, orient, usedL, pos = normalize_to_reference(
            seq, ref_full, args.anchor_len, args.min_anchor, args.step
        )

        out_fp = out_dir / fp.name
        if norm is None:
            # Fallback: write the original sequence and mark as not found
            write_fasta(out_fp, sample, seq)
            report.append(
                f"{fp.name}\t{sample}\tanchor_not_found\t{orient}\t{usedL}\t{pos}\t{len(seq)}\t{out_fp.name}"
            )
        else:
            write_fasta(out_fp, sample, norm)
            report.append(
                f"{fp.name}\t{sample}\tok\t{orient}\t{usedL}\t{pos}\t{len(norm)}\t{out_fp.name}"
            )

    (out_dir / "normalization_report.tsv").write_text("\n".join(report) + "\n")


if __name__ == "__main__":
    main()
