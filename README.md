# Restriction Enzyme Mapper

A Python command-line tool that finds restriction enzyme recognition sites in a DNA sequence, calculates where each enzyme cuts, predicts the fragment sizes of a digest, and draws a site map and a virtual agarose gel.

Built by **Garv Gupta** (B.Tech Biotechnology)  Its cut positions are cross-checked against Biopython's `Bio.Restriction`.



*Example: EcoRV + PstI + BglII digest of the 884 bp sequence in `data/sample.fasta`.*

## Features

- Reads standard FASTA files (multi-line sequences, lowercase, multiple records)
- Supports **IUPAC degenerate bases** in recognition sites (for example `GTYRAC` for HincII)
- Finds **overlapping** sites
- Reports the recognition window and the exact cut position for each hit
- Calculates fragment sizes for **single or combined (multi-enzyme) digests**
- Handles **linear and circular** DNA, including the wrap-around fragment of a plasmid
- Exports a CSV report
- Optional site map and virtual gel (PNG)
- Unit tests, plus validation against Biopython

## Background

Type II restriction enzymes recognise short DNA sequences (usually 4 to 8 bp, often palindromic) and cut at a fixed position inside or next to them. For example, EcoRI recognises `GAATTC` and cuts between the G and the first A (`G^AATTC`). Some enzymes tolerate variable bases at certain positions, which is written with IUPAC codes. HincII recognises `GTYRAC`, where Y is C or T and R is A or G. Cutting a sequence with one or more enzymes produces fragments whose sizes can be predicted, which is how restriction digests and agarose gels are planned in cloning work.

## Installation

```
git clone https://github.com/<your-username>/restriction-mapper.git
cd restriction-mapper
pip install -r requirements.txt
```

Python 3.9 or newer is recommended. `matplotlib` is needed for the plots, `biopython` for the validation script, and `pytest` for the tests.

## Usage

Run with the defaults (all enzymes in the database, linear DNA, sample sequence):

```
python restriction_mapper.py
```

Choose enzymes and draw the plots:

```
python restriction_mapper.py --enzymes EcoRV PstI BglII --plots
```

Treat the DNA as a circular plasmid:

```
python restriction_mapper.py --enzymes EcoRV PstI --topology circular
```

Use your own FASTA file and output location:

```
python restriction_mapper.py -i my_plasmid.fasta -o results/my_report.csv --plots
```

### Options

| Option | Description |
|---|---|
| `-i`, `--input` | Input FASTA file (default: `data/sample.fasta`) |
| `--db` | Enzyme database in JSON (default: `enzymes_db.json`) |
| `--enzymes` | Enzymes to use together in one digest (default: all in the database) |
| `--topology` | `linear` (default) or `circular` |
| `-o`, `--output` | Output CSV file (default: `results/digestion_report.csv`) |
| `--plots` | Also save a site map and a virtual gel as PNG files next to the CSV |

On Windows, use `py` instead of `python`.

### Example output

```
JS684847.1  length = 884 bp  (linear)

Enzyme   Site      Start   End  Cut after
BglII    AGATCT       27    32         27
EcoRV    GATATC      272   277        274
PstI     CTGCAG      752   757        756
EcoRV    GATATC      853   858        855

Fragments (5): [482, 247, 99, 29, 27]  sum = 884
```

The same EcoRV + PstI digest on circular DNA gives 3 fragments, `[482, 303, 99]`, because the two ends of the sequence join into one fragment.

`Cut after` means the top strand is cut after that base (1-based).

### CSV report

```
sequence,enzyme,site,start,end,cut_after
JS684847.1,BglII,AGATCT,27,32,27
JS684847.1,EcoRV,GATATC,272,277,274
JS684847.1,PstI,CTGCAG,752,757,756
JS684847.1,EcoRV,GATATC,853,858,855
```

## Enzyme database

Enzymes are stored in `enzymes_db.json`:

```json
"EcoRI": {"site": "GAATTC", "cut_top": 1, "cut_bottom": 5}
```

- `site`: recognition sequence (IUPAC codes allowed)
- `cut_top`: number of bases from the start of the site to the top-strand cut
- `cut_bottom`: the equivalent position on the bottom strand, kept for reference (fragment sizes use the top-strand cut)

To add an enzyme, add a line with its site and cut offsets. The database currently holds 14 enzymes: EcoRI, BamHI, HindIII, XbaI, BglII, SalI, NcoI, PstI, KpnI, EcoRV, SmaI, HaeIII, HincII and AvaII.

## Testing and validation

Run the unit tests:

```
python -m pytest
```

They cover IUPAC translation, overlapping sites, cut positions, linear and circular fragment sizes (including the check that fragments always add up to the sequence length), and FASTA parsing.

Cross-check against Biopython:

```
python validate_biopython.py
```

This compares the cut positions found by this tool with `Bio.Restriction` for every enzyme in the database and prints `OK` or `MISMATCH` for each.

## Project structure

```
restriction-mapper/
├── restriction_mapper.py        # main command-line tool
├── enzymes_db.json              # enzyme recognition sites and cut offsets
├── test_restriction_mapper.py   # unit tests
├── validate_biopython.py        # comparison with Biopython
├── requirements.txt
├── LICENSE
├── data/
│   └── sample.fasta             # test sequence (GenBank JS684847.1)
└── results/                     # example CSV report, map and gel
```

## Limitations

- Only the forward strand is searched. This is correct for palindromic sites, but non-palindromic enzymes (for example BsaI or BsmBI) would also need a reverse-strand search.
- The database has 14 common enzymes, not the full REBASE collection.
- Methylation sensitivity (Dam/Dcm) and star activity are not modelled.
- Only a linear site map is drawn, with no circular plasmid map.
- Cut positions follow the top strand. Overhang lengths are stored in the database but not reported.

## Roadmap

- Reverse-strand scanning for non-palindromic enzymes
- Circular plasmid map
- Loading enzymes directly from a REBASE export
- Reporting overhang type (5', 3' or blunt)
- Multiple digests side by side in one gel image

## Data

`data/sample.fasta` is the public GenBank record JS684847.1 (*Streblospio benedicti* genomic survey sequence from a pUC19 library clone, 884 bp), used as a test input.

## Author

Garv Gupta, B.Tech Biotechnology
