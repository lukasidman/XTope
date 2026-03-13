#!/usr/bin/env python3
"""
Generate a large benchmark antigen dataset for timing comparisons
between different alignment methods.

Produces ~5000 sequences in the same semicolon-delimited format as
benchmark_antigens.csv, with planted homologous families scaled up
to ensure cross-reactive pairs exist at realistic density.

Usage:
    python generate_large_benchmark.py [--size 5000] [--output large_benchmark.csv]
"""

import argparse
import csv
import random
import sys

random.seed(2026)

AA = "ACDEFGHIKLMNPQRSTVWY"
TAG = "MHHHHHHGSSG"


def random_seq(length: int) -> str:
    return "".join(random.choices(AA, k=length))


def mutate(seq: str, rate: float) -> str:
    """Point-mutate each position independently at the given rate."""
    return "".join(
        c if random.random() > rate else random.choice(AA)
        for c in seq
    )


def insert_motif(seq: str, motif: str, position: int) -> str:
    """Splice a motif into a sequence at the given position."""
    return seq[:position] + motif + seq[position + len(motif):]


def make_family_fulllength(name: str, n_members: int, ancestor_len: int,
                           min_rate: float = 0.03, max_rate: float = 0.18) -> list[dict]:
    """Full-length homology family — all members derived from one ancestor."""
    ancestor = random_seq(ancestor_len)
    records = []
    for i in range(n_members):
        rate = min_rate + (max_rate - min_rate) * i / max(n_members - 1, 1)
        variant = mutate(ancestor, rate)
        records.append({
            "name": f"{name}_{i+1:02d}",
            "sequence": TAG + variant,
        })
    return records


def make_family_motif(name: str, n_members: int, motif: str,
                      flank_range: tuple[int, int] = (30, 70),
                      motif_mut_rate: float = 0.08) -> list[dict]:
    """Local-motif family — shared conserved core in random flanking sequence."""
    records = []
    for i in range(n_members):
        left_len = random.randint(*flank_range)
        right_len = random.randint(*flank_range)
        left = random_seq(left_len)
        right = random_seq(right_len)
        core = mutate(motif, motif_mut_rate) if i > 0 else motif
        seq = left + core + right
        records.append({
            "name": f"{name}_{i+1:02d}",
            "sequence": TAG + seq,
        })
    return records


def make_family_epitope(name: str, n_members: int, epitope: str,
                        bg_range: tuple[int, int] = (55, 90)) -> list[dict]:
    """Short-epitope family — small shared peptide in diverse backgrounds."""
    records = []
    for i in range(n_members):
        bg_len = random.randint(*bg_range)
        bg = random_seq(bg_len)
        pos = random.randint(10, bg_len - len(epitope) - 10)
        seq = insert_motif(bg, epitope, pos)
        records.append({
            "name": f"{name}_{i+1:02d}",
            "sequence": TAG + seq,
        })
    return records


def make_family_convergent(name: str, n_members: int, length: int = 80) -> list[dict]:
    """
    Physicochemical-convergent family — different sequences, similar
    hydrophobicity/charge pattern. Uses alternating hydrophobic/polar blocks.
    """
    hydrophobic = "AVILMF"
    polar = "DEKRNQ"
    pattern_block = 8  # residues per block

    records = []
    for i in range(n_members):
        seq = []
        for b in range(length // pattern_block + 1):
            pool = hydrophobic if b % 2 == 0 else polar
            seq.extend(random.choices(pool, k=pattern_block))
        records.append({
            "name": f"{name}_{i+1:02d}",
            "sequence": TAG + "".join(seq[:length]),
        })
    return records


def make_family_dual_domain(name: str, n_members: int,
                            domain_a: str, domain_b: str,
                            linker_range: tuple[int, int] = (15, 35),
                            mut_rate: float = 0.10) -> list[dict]:
    """Two conserved domains separated by a variable linker."""
    records = []
    for i in range(n_members):
        da = mutate(domain_a, mut_rate) if i > 0 else domain_a
        db = mutate(domain_b, mut_rate) if i > 0 else domain_b
        linker_len = random.randint(*linker_range)
        left_flank = random_seq(random.randint(10, 25))
        right_flank = random_seq(random.randint(10, 25))
        seq = left_flank + da + random_seq(linker_len) + db + right_flank
        records.append({
            "name": f"{name}_{i+1:02d}",
            "sequence": TAG + seq,
        })
    return records


def generate_dataset(total: int) -> list[dict]:
    """Build the full dataset: planted families + random background."""
    planted: list[dict] = []

    # ── Scaled-up versions of the 8 original benchmark families ──

    # 1. Full-length homology (multiple sizes)
    for j in range(4):
        ancestor_len = random.randint(60, 120)
        n = random.randint(6, 12)
        planted.extend(make_family_fulllength(
            f"CYTK{j+1:02d}", n, ancestor_len, 0.03, 0.20
        ))

    # 2. Local motif families
    motifs = [
        "PENILFLDAHYCLSAGELKWLTKCPD",
        "RVWQTDMCYSHGPNFKEL",
        "DCFGKLAPENVHWTSYMRI",
        "AQWSNCGDTIRFHPKELVM",
        "GHTMWKRSDCEPFVNLIYA",
    ]
    for j, motif in enumerate(motifs):
        n = random.randint(5, 10)
        planted.extend(make_family_motif(f"KIN{j+1:02d}", n, motif))

    # 3. Short epitope families
    epitopes = [
        "RPCHQFNV",
        "DCWGKYME",
        "FHTPSNVL",
        "GKWRCEDY",
        "NMPQHFAT",
        "VCSRKDWI",
    ]
    for j, ep in enumerate(epitopes):
        n = random.randint(4, 8)
        planted.extend(make_family_epitope(f"EPI{j+1:02d}", n, ep))

    # 4. Convergent families (physicochemical only)
    for j in range(3):
        length = random.randint(60, 100)
        n = random.randint(4, 7)
        planted.extend(make_family_convergent(f"CONV{j+1:02d}", n, length))

    # 5. Dual-domain families
    domains = [
        ("CVWTHKGDNECRY", "FPAQMLSVNHKWT"),
        ("EDCRFGHKNPWVY", "AGILMSTCDEQRK"),
        ("WYHFCKNDRSEPM", "TILVAGSQNHKDE"),
    ]
    for j, (da, db) in enumerate(domains):
        n = random.randint(5, 8)
        planted.extend(make_family_dual_domain(f"DUAL{j+1:02d}", n, da, db))

    # 6. Divergence gradient families — single ancestor, wide mutation range
    for j in range(3):
        ancestor_len = random.randint(70, 110)
        planted.extend(make_family_fulllength(
            f"GRAD{j+1:02d}", 8, ancestor_len, 0.02, 0.55
        ))

    # ── Fill remaining slots with random background sequences ──
    n_bg = total - len(planted)
    if n_bg < 0:
        print(f"Warning: planted families ({len(planted)}) exceed requested "
              f"total ({total}). Producing {len(planted)} sequences.", file=sys.stderr)
        n_bg = 0

    background = []
    for i in range(n_bg):
        length = random.randint(40, 140)
        background.append({
            "name": f"BG_{i+1:04d}",
            "sequence": TAG + random_seq(length),
        })

    all_records = planted + background
    random.shuffle(all_records)
    return all_records


def main():
    parser = argparse.ArgumentParser(
        description="Generate a large benchmark antigen dataset"
    )
    parser.add_argument(
        "--size", type=int, default=5000,
        help="Total number of sequences (default: 5000)"
    )
    parser.add_argument(
        "--output", type=str, default="large_benchmark.csv",
        help="Output CSV filename (default: large_benchmark.csv)"
    )
    args = parser.parse_args()

    records = generate_dataset(args.size)

    with open(args.output, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter=";")
        writer.writerow(["", "name", "sequence"])  # match benchmark format
        for rec in records:
            writer.writerow(["", rec["name"], rec["sequence"]])

    # Summary
    planted = [r for r in records if not r["name"].startswith("BG_")]
    families: dict[str, int] = {}
    for r in planted:
        fam = r["name"].rsplit("_", 1)[0]
        families[fam] = families.get(fam, 0) + 1

    print(f"Generated {len(records)} sequences → {args.output}")
    print(f"  Planted: {len(planted)} sequences in {len(families)} families")
    print(f"  Background: {len(records) - len(planted)} random sequences")
    print(f"\nFamily breakdown:")
    for fam, count in sorted(families.items()):
        print(f"  {fam}: {count} members")


if __name__ == "__main__":
    main()
