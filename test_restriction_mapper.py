"""Unit tests. Run with:  py -m pytest

Author: Garv Gupta, B.Tech Biotechnology
"""
from restriction_mapper import (
    BASE_DIR, read_fasta, load_enzymes, motif_to_regex,
    find_sites, find_cuts, fragment_sizes, digest,
)

ENZYMES = load_enzymes(BASE_DIR / "enzymes_db.json")
SAMPLE = read_fasta(BASE_DIR / "data" / "sample.fasta")["JS684847.1"]


# ---- IUPAC / motif matching ----

def test_iupac_regex_translation():
    assert motif_to_regex("GTYRAC") == "GT[CT][AG]AC"
    assert motif_to_regex("GGWCC") == "GG[AT]CC"


def test_degenerate_motif_matches_all_variants():
    # GTYRAC should match GTCAAC, GTCGAC, GTTAAC, GTTGAC
    for variant in ("GTCAAC", "GTCGAC", "GTTAAC", "GTTGAC"):
        assert find_sites(variant, "GTYRAC") == [0]
    assert find_sites("GTAAAC", "GTYRAC") == []     # A is not a pyrimidine


def test_overlapping_sites_are_found():
    assert find_sites("GGGG", "GGG") == [0, 1]


# ---- cut positions ----

def test_ecori_cut_position():
    cuts = find_cuts("AAAGAATTCAAA", {"EcoRI": ENZYMES["EcoRI"]})
    assert len(cuts) == 1
    assert cuts[0]["start"] == 4 and cuts[0]["end"] == 9
    assert cuts[0]["cut_after"] == 4            # G^AATTC: cut after the G


def test_no_sites_gives_no_cuts():
    assert find_cuts("AAAAAAAAAA", ENZYMES) == []


# ---- fragment sizes ----

def test_linear_fragments_sum_to_length():
    sizes = fragment_sizes([100, 400], 1000)
    assert sorted(sizes) == [100, 300, 600]
    assert sum(sizes) == 1000


def test_circular_wraparound_fragment():
    sizes = fragment_sizes([100, 400], 1000, circular=True)
    assert sorted(sizes) == [300, 700]          # n cuts -> n fragments
    assert sum(sizes) == 1000


def test_circular_single_cut_linearizes():
    assert fragment_sizes([250], 1000, circular=True) == [1000]


def test_no_cuts_returns_whole_sequence():
    assert fragment_sizes([], 1000) == [1000]
    assert fragment_sizes([], 1000, circular=True) == [1000]


# ---- full digests on the sample sequence ----

def test_sample_length():
    assert len(SAMPLE) == 884


def test_ecorv_digest_on_sample():
    _, sizes = digest(SAMPLE, ENZYMES, ["EcoRV"])
    assert sizes == [581, 274, 29]


def test_double_digest_on_sample():
    _, sizes = digest(SAMPLE, ENZYMES, ["EcoRV", "PstI"])
    assert sizes == [482, 274, 99, 29]


def test_fragments_always_sum_to_length():
    for names in (["EcoRV"], ["EcoRV", "PstI"], ["EcoRV", "PstI", "BglII"], list(ENZYMES)):
        for circular in (False, True):
            _, sizes = digest(SAMPLE, ENZYMES, names, circular)
            assert sum(sizes) == len(SAMPLE)


# ---- FASTA parsing ----

def test_read_fasta_multiline_and_lowercase(tmp_path):
    fasta = tmp_path / "t.fasta"
    fasta.write_text(">seq1 description\nacgt\nACGT\n>seq2\nTTTT\n")
    records = read_fasta(fasta)
    assert records == {"seq1": "ACGTACGT", "seq2": "TTTT"}
