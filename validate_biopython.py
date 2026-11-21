"""
Cross-check this tool against Biopython's Bio.Restriction.

For every enzyme in the database, compare the cut positions found by
restriction_mapper.py with those reported by Biopython on the sample
sequence. Requires:  pip install biopython

Author: Garv Gupta, B.Tech Biotechnology
"""
from Bio import Restriction
from Bio.Seq import Seq

from restriction_mapper import BASE_DIR, read_fasta, load_enzymes, find_cuts

records = read_fasta(BASE_DIR / "data" / "sample.fasta")
enzymes = load_enzymes(BASE_DIR / "enzymes_db.json")

all_ok = True
for rec_name, seq in records.items():
    print(f"Sequence: {rec_name} ({len(seq)} bp)\n")
    print(f"{'Enzyme':<9}{'This tool':<22}{'Biopython':<22}Result")
    for enz_name, info in enzymes.items():
        mine = sorted(c["cut_after"] for c in find_cuts(seq, {enz_name: info}))
        bio_enzyme = getattr(Restriction, enz_name)
        # Biopython reports the 1-based position of the first base AFTER the cut,
        # so subtract 1 to get "cut after base N", the convention used here
        theirs = sorted(p - 1 for p in bio_enzyme.search(Seq(seq), linear=True))
        ok = mine == theirs
        all_ok &= ok
        print(f"{enz_name:<9}{str(mine):<22}{str(theirs):<22}{'OK' if ok else 'MISMATCH'}")

print("\nAll enzymes agree with Biopython." if all_ok else
      "\nSome enzymes differ - check the rows marked MISMATCH.")
