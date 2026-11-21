"""
Restriction Enzyme Mapper
Find restriction enzyme recognition sites in a DNA sequence (FASTA),
calculate cut positions and digest fragment sizes, and export a CSV
report plus an optional site map and virtual gel.

Author: Garv Gupta, B.Tech Biotechnology
"""
import re
import csv
import json
import argparse
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent

# IUPAC nucleotide codes -> regular expression character classes
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}


def read_fasta(path):
    """Read a FASTA file into {record_id: UPPERCASE sequence}."""
    records, name, seq = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    records[name] = "".join(seq).upper()
                name, seq = line[1:].split()[0], []
            elif line:
                seq.append(line)
    if name:
        records[name] = "".join(seq).upper()
    return records


def load_enzymes(path):
    with open(path) as f:
        return json.load(f)


def motif_to_regex(motif):
    """Translate an IUPAC recognition motif into a regex pattern."""
    return "".join(IUPAC[base] for base in motif)


def find_sites(seq, motif):
    """0-based start positions of every match, including overlapping ones."""
    pattern = re.compile(f"(?=({motif_to_regex(motif)}))")
    return [m.start() for m in pattern.finditer(seq)]


def find_cuts(seq, enzymes):
    """Find every recognition site and the top-strand cut position.

    'cut_after' = the top strand is cut after this base (1-based).
    """
    results = []
    for name, info in enzymes.items():
        site = info["site"]
        for start in find_sites(seq, site):
            results.append({
                "enzyme": name,
                "site": site,
                "start": start + 1,
                "end": start + len(site),
                "cut_after": start + info["cut_top"],
            })
    return sorted(results, key=lambda r: r["start"])


def fragment_sizes(cut_positions, seq_len, circular=False):
    """Fragment sizes (bp) produced by cutting at the given positions."""
    cuts = sorted(set(cut_positions))
    if not cuts:
        return [seq_len]
    if circular:
        sizes = [b - a for a, b in zip(cuts, cuts[1:])]
        sizes.append(seq_len - cuts[-1] + cuts[0])      # wrap-around fragment
    else:
        points = [0] + cuts + [seq_len]
        sizes = [b - a for a, b in zip(points, points[1:])]
    return sizes


def digest(seq, enzymes, names, circular=False):
    """Digest a sequence with one or more enzymes at the same time."""
    chosen = {n: enzymes[n] for n in names}
    cuts = find_cuts(seq, chosen)
    positions = [c["cut_after"] for c in cuts]
    sizes = sorted(fragment_sizes(positions, len(seq), circular), reverse=True)
    return cuts, sizes


def write_csv(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["sequence", "enzyme", "site", "start", "end", "cut_after"]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def plot_map(name, seq_len, cuts, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    enzyme_names = sorted({c["enzyme"] for c in cuts})
    colors = {n: plt.cm.tab10(i % 10) for i, n in enumerate(enzyme_names)}

    fig, ax = plt.subplots(figsize=(11, 3.2))
    ax.hlines(0, 0, seq_len, color="black", linewidth=5)
    for i, c in enumerate(cuts):
        x = c["cut_after"]
        h = 0.6 + 0.4 * (i % 3)          # stagger heights so labels don't overlap
        col = colors[c["enzyme"]]
        ax.vlines(x, 0, h, color=col, linewidth=1.5)
        ax.text(x, h + 0.05, f"{c['enzyme']}\n{x}", ha="center",
                va="bottom", fontsize=8, color=col)

    ax.set_xlim(-seq_len * 0.03, seq_len * 1.03)
    ax.set_ylim(-0.4, 2.4)
    ax.set_yticks([])
    ax.set_xlabel("Position (bp)")
    ax.set_title(f"{name} restriction map ({seq_len} bp, {len(cuts)} sites)")
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def plot_gel(name, sizes, path, label):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ladder = [1500, 1200, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100]

    fig, ax = plt.subplots(figsize=(4, 6))
    ax.set_facecolor("black")
    for s in ladder:
        ax.hlines(s, 0.6, 1.4, color="white", linewidth=3)
    for s in sizes:
        ax.hlines(s, 1.6, 2.4, color="#ffd54f", linewidth=4)
        ax.text(2.5, s, f"{s}", color="white", va="center", fontsize=8)

    ax.set_yscale("log")
    ax.set_yticks(ladder)
    ax.set_yticklabels([str(s) for s in ladder])
    ax.minorticks_off()
    ax.set_ylim(min(sizes + [100]) * 0.7, max(sizes + [1500]) * 1.2)
    ax.set_xlim(0, 3.2)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["Ladder", label], fontsize=8)
    ax.set_ylabel("Size (bp)")
    ax.set_title(f"{name}\nvirtual gel", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def parse_args():
    p = argparse.ArgumentParser(description="Restriction enzyme mapper")
    p.add_argument("-i", "--input", default=str(BASE_DIR / "data" / "sample.fasta"),
                   help="input FASTA file")
    p.add_argument("--db", default=str(BASE_DIR / "enzymes_db.json"),
                   help="enzyme database (JSON)")
    p.add_argument("--enzymes", nargs="+",
                   help="enzyme names to use (default: all in the database)")
    p.add_argument("--topology", choices=["linear", "circular"], default="linear")
    p.add_argument("-o", "--output", default=str(BASE_DIR / "results" / "digestion_report.csv"),
                   help="output CSV file")
    p.add_argument("--plots", action="store_true",
                   help="also save a site map and a virtual gel (PNG) next to the CSV")
    return p.parse_args()


def main():
    args = parse_args()
    records = read_fasta(args.input)
    enzymes = load_enzymes(args.db)

    names = args.enzymes or list(enzymes)
    unknown = [n for n in names if n not in enzymes]
    if unknown:
        raise SystemExit(f"Unknown enzyme(s): {', '.join(unknown)}")

    circular = args.topology == "circular"
    out_dir = Path(args.output).parent
    rows = []

    for name, seq in records.items():
        cuts, sizes = digest(seq, enzymes, names, circular)
        print(f"{name}  length = {len(seq)} bp  ({args.topology})\n")
        print(f"{'Enzyme':<9}{'Site':<9}{'Start':>6}{'End':>6}{'Cut after':>11}")
        for r in cuts:
            print(f"{r['enzyme']:<9}{r['site']:<9}{r['start']:>6}{r['end']:>6}{r['cut_after']:>11}")
            rows.append({"sequence": name, **r})
        print(f"\nFragments ({len(sizes)}): {sizes}  sum = {sum(sizes)}")

        if args.plots:
            out_dir.mkdir(parents=True, exist_ok=True)
            map_path = out_dir / f"{name}_map.png"
            gel_path = out_dir / f"{name}_gel.png"
            plot_map(name, len(seq), cuts, map_path)
            plot_gel(name, sizes, gel_path, " + ".join(names) if args.enzymes else "all")
            print(f"Plots saved: {map_path.name}, {gel_path.name}")

    write_csv(args.output, rows)
    print(f"\nCSV saved to {args.output}")


if __name__ == "__main__":
    main()
