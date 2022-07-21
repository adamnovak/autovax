#!/usr/bin/env python3.9

import logging
import random
import re
import sys
from argparse import Namespace

from Bio import SeqIO
from Bio.Seq import Seq, transcribe
from Bio.SeqRecord import SeqRecord
from Bio.pairwise2 import align, format_alignment
from Bio.Align import MultipleSeqAlignment

# Codon Harmony now seems to live at https://github.com/ArsenicMeatball/codon-harmony
# See <https://github.com/jocelynnpearl/codon_tools/blob/4b84b851a14202724f5c3605b0878d77c5d598d1/codon_harmony/codon_harmony.py>
# To use Codon Harmony, we need an old BioPython that still has Alphabet.
# See <https://biopython.org/wiki/Alphabet>
from codon_harmony.util import codon_use
from codon_harmony.data import GC_content, RibosomeBindingSites, RestrictionEnzymes
from codon_harmony.codon_harmony import _harmonize_sequence as harmonize_sequence

logger = logging.getLogger(__name__)

def align_sequence_to_sequence(needle: Seq, haystack: Seq) -> MultipleSeqAlignment:
    """
    Align the neelde sequence to the haystack sequence, with free end gaps.
    """
    
    # Find the optimal global alignment with free gaps in the needle sequence
    alignments = align.globalms(
        needle,
        haystack,
        match=6,
        mismatch=-1,
        open=-4,
        extend=-1,
        penalize_end_gaps=(False, True),
        one_alignment_only=True
    )

    if len(alignments) == 0:
        raise RuntimeError(f"Cannot find sequence")
    
    alignment = alignments[0]
    return alignment

def find_sequence_in_sequence(needle: Seq, haystack: Seq, max_length_error: float = 0.1) -> tuple[tuple[int, int], list[tuple[int, int, int]]]:
    """
    Find the [start, end) occupied by the needle sequence in the
    haystack sequence. Also returns regions of correspondence between the
    needle and the identified range of the haystack.
    
    Throws an error if the length difference is more than max_length_error of
    the needle's length.
    """
    
    alignment = align_sequence_to_sequence(needle, haystack)
    gapped_needle = alignment[0]
    gapped_haystack = alignment[1]
    
    # Count the leading gap characters
    leading_gap = 0
    while leading_gap < len(gapped_needle) and gapped_needle[leading_gap] == '-':
        leading_gap += 1
        
    # Count the trailing gap characters
    trailing_gap = 0
    while trailing_gap < len(gapped_needle) and gapped_needle[-(trailing_gap + 1)] == '-':
        trailing_gap += 1
        
    # The range occupied by the needle is all of the haystack except the parts
    # skipped by the leading and trailing gaps   
    overlapped_length = len(haystack) - leading_gap - trailing_gap
    
    length_change = abs(overlapped_length - len(needle))
    if length_change/len(needle) > max_length_error:
        # This is probably not the right thing.
        raise RuntimeError(f"Length difference of {length_change} bp for {len(needle)} bp query is excessive")
        
    # To translate from haystack coordinates relative to the start into needle coordinates, we need to know the matched ranges.
    # Each matched range is a haystack start offset, a needle start offset, and a length
    matched_ranges = []
    current_range = None
        
    needle_cursor = 0
    haystack_cursor = 0
    
    for i in range(leading_gap, len(gapped_needle) - trailing_gap):
        if gapped_needle[i] == '-':
            # Gap
            if current_range is not None:
                matched_ranges.append(tuple(current_range))
                current_range = None
            # Only move in haystack
            haystack_cursor += 1
        elif gapped_haystack[i] == '-':
            # Gap
            if current_range is not None:
                matched_ranges.append(tuple(current_range))
                current_range = None
            # Only move in needle
            needle_cursor += 1
        else:
            # Match/mismatch, move in both
            if current_range is None:
                current_range = [needle_cursor, haystack_cursor, 1]
            else:
                current_range[2] += 1
        
            # Move in both
            needle_cursor += 1
            haystack_cursor += 1
            
    if current_range is not None:
        matched_ranges.append(tuple(current_range))
        current_range = None
        
    # And the needle starts at the end of the part of the haystack skipped by
    # the leading gap.
    return (leading_gap, leading_gap + overlapped_length), matched_ranges
    
def extract_spike(spike_dna: SeqRecord, target_dna: SeqRecord) -> SeqRecord:
    """
    Find and extract the part of a full SARS-Cov-2 sequence corresponding to
    the given spike ORF template.
    """
    
    bounds, _ = find_sequence_in_sequence(spike_dna.seq, target_dna.seq)
    found_dna = target_dna[bounds[0]:bounds[1]]
    return found_dna
    
def dna_to_protein(template_dna: SeqRecord) -> Seq:
    """
    Convert DNA to a protein. Raises an error if the protein doesn't look plausible.
    """
    
    protein_seq = transcribe(template_dna.seq).translate()
    if len(protein_seq) * 3 != len(template_dna.seq) or protein_seq[0] != 'M' or protein_seq[-1] != '*' or '*' in protein_seq[:-1]:
        raise RuntimeError(f"Cannot operate on broken protein: {protein_seq}")
    return protein_seq
    
def dna_regions_to_protein_regions(dna_matches: list[tuple[int, int, int]]) -> list[tuple[int, int, int]]:
    """
    Change matched regions in DNA space into matched regions in protein space.
    
    Raises an error if there are frame shifts.
    """
    
    def change_number(number: int) -> int:
        if number % 3 != 0:
            raise RuntimeError("Frame-shift detected")
        return number / 3
    
    return [tuple([change_number(x) for x in t]) for t in dna_matches]
    
def lift_coordinate(coordinate: int, matches: list[tuple[int, int, int]]) -> int:
    """
    Given a coordinate in the canonical sequence, and a list of tuples of
    (canonical offset, found offset, length) matches, lift over the coordinate
    from the canonical sequence to the found sequence.
    
    Raises an error if it cannot be lifted over.
    """
    
    # TODO: Binary search!
    for source_start, dest_start, length in matches:
        if source_start <= coordinate and source_start + length > coordinate:
            return dest_start + (coordinate - source_start)
    
    raise RuntimeError(f"Coordinate {coordinate} does not fall in a match block in {matches}")

def remove_furin_cleavage(protein_seq: Seq) -> Seq:
    """
    Remove the furin cleavage site (RRAR) and replace it with GSAS as
    recommended in <https://www.nature.com/articles/s41467-022-30878-4>.
    """
    
    FURIN_SITE = "RRAR"
    REPLACEMENT = "GSAS"
    
    if protein_seq.count(FURIN_SITE) != 1:
        raise RuntimeError("No unique furin cleavage site found to replace")
    
    # Make it mutable
    new_seq = protein_seq.tomutable()
    
    # Find and overwrite the furin site
    index = protein_seq.find(FURIN_SITE)
    for i in range(len(FURIN_SITE)):
        new_seq[index + i] = REPLACEMENT[i]
    
    return new_seq.toseq()
    
def hexapro(protein_seq: Seq, protein_matches: list[tuple[int, int, int]]) -> Seq:
    """
    Replace six AAs with Proline as recommended in
    <https://www.science.org/doi/10.1126/science.abd0826>. Actual substitution locations from <https://www.addgene.org/154754/>.
    """
    
    S2P_SUBSTITUTIONS = ["p.K986P", "p.V987P"] 
    HEXAPRO_SUBSTITUTIONS = ["p.F817P", "p.A892P", "p.A899P", "p.A942P"]
    
    # We would parse and apply these with
    # <https://hgvs.readthedocs.io/en/stable/installation.html> but that
    # apparently needs Postgres headers to... parse strings. So we do it
    # ourselves.
    
    new_seq = protein_seq.tomutable()
    
    for hgvs_spec in S2P_SUBSTITUTIONS + HEXAPRO_SUBSTITUTIONS:
        match = re.match("p\\.([A-Z])([0-9]+)([A-Z])", hgvs_spec)
        from_aa = match.group(1)
        to_aa = match.group(3)
        one_based_position = int(match.group(2))
        
        # Convert position to this sequence's coordinates from the alignment
        lifted_position = lift_coordinate(one_based_position - 1, protein_matches)
        
        if new_seq[lifted_position] != from_aa:
            # We don't have the expected AA here. Seems suspicious.
            logger.warning("Applying %s at %s where we actually have %s", hgvs_spec, lifted_position, new_seq[lifted_position])
        # Use the lifted-over position to rewrite the appropriate amino acid
        new_seq[lifted_position] = to_aa
        
    return new_seq.toseq()
    
def protein_to_optimized_dna(protein_seq: Seq, mode: str = "normal") -> Seq:
    """
    Translate a protein sequence back to codon-optimized DNA for humans.
    
    The sequence must end in a stop codon (*)
    """
    
    if protein_seq[-1] != '*':
        raise RuntimeError("Cannot back-translate unterminated protien")
    
    HUMAN_TAXID = "9606"
    RESTRICTION_ENZYMES = ["NdeI", "XhoI", "HpaI", "PstI", "EcoRV", "NcoI", "BamHI"]
    STOP_CODON = "TAG"
    
    # Duplicate logic from codon_harmony because we don't want file IO.
    codon_use_table, host_profile, codon_relative_adativeness = codon_use.host_codon_usage(HUMAN_TAXID)
    rest_enz = RestrictionEnzymes(RESTRICTION_ENZYMES)
    fake_args = Namespace()
    # Default cycles and inner_cycles are 10 each but that takes forever
    fake_args.cycles = 3 if mode == "normal" else 1
    fake_args.inner_cycles = 3 if mode == "normal" else 1
    fake_args.max_relax = 0.1
    fake_args.start_sites = False
    fake_args.local_homopolymer_threshold = 4
    fake_args.splice_sites = True
    # codon_harmony can't handle stop codons
    result = harmonize_sequence(
        SeqRecord(protein_seq[:-1]),
        fake_args,
        codon_use_table,
        host_profile,
        codon_relative_adativeness,
        rest_enz,
        0
    )
    if 'dna' not in result:
        raise RuntimeError("Could not create codon sequence")
    return Seq(result['dna'].split('\n', 1)[1].replace('\n', '') + STOP_CODON)
    
def main():

    if len(sys.argv) != 2:
        print("Autovax: Automatically design a SARS-CoV-2 vaccine candidate from a sequenced sample.")
        print()
        print("Attempts to follow the methods of <https://www.nature.com/articles/s41467-022-30878-4>, using commercial synthesis. Produces instructions for ordering mRNA and encapsulating into nanoparticles.")
        print()
        print(f"Usage: {sys.argv[0]} <sampled-sars-cov-2-genome.fa>")
        print()
        sys.exit(1)

    # Set to "fast" or "normal".
    # "fast" will skip some computation to exercise the code faster.
    mode = "normal"

    logging.basicConfig(level=logging.INFO)
    # TODO: Suppress Codon Harmony's logger.
    
    # Be deterministic
    random.seed(1)
    
    logger.info("Load sample sequence...")
    target_dna = SeqIO.read(sys.argv[1], 'fasta')

    logger.info("Load canonical spike DNA...")
    # This was extracted from the UCSC genome browser
    spike_dna = SeqIO.read('spike.fa', 'fasta')
    logger.info("Translate canonical spike to protein...")
    canonical_protein_seq = dna_to_protein(spike_dna)
    
    logger.info("Extract sample's spike DNA...")
    if mode == "normal":
        # Find the spike sequence
        found_dna = extract_spike(spike_dna, target_dna)
    else:
        # Fake it
        found_dna = spike_dna
        dna_matches = [(0, 0, len(spike_dna.seq))]
        
    # Apply transformations
    logger.info("Translate sample spike to protein...")
    protein_seq = dna_to_protein(found_dna)
    
    logger.info("Establish canonical coordinate system on sample spike protein...")
    _, protein_matches = find_sequence_in_sequence(canonical_protein_seq, protein_seq)
    
    logger.info("Remove Furin cleavage site from sample spike protein...")
    protein_seq = remove_furin_cleavage(protein_seq)
    logger.info("Apply Hexapro Prolines to sample spike protein...")
    protein_seq = hexapro(protein_seq, protein_matches)
    
    logger.info("Back-translate stabilized sample spike protein to codon-optimized DNA...")
    dna_seq = protein_to_optimized_dna(protein_seq, mode=mode)
    
    logger.info("Prepare RNA order...")
    
    print()
    print()
    print("RNA DETAILS")
    print()
    print(f"RNA Name: Vaccine candidate derived from {found_dna.id}") 
    print("Modifications: 'N1-Methylpseudouridin' [sic]")
    print("Type of RNA: Capped mRNA")
    print("Template Source: TriLink Supplied")
    print(f"ORF Sequence: {dna_seq}")
    print()
    print("MANUFACTURING PROCEDURE")
    print()
    print("1. Enter above details into <https://www.trilinkbiotech.com/mrna.html> and order.")
    print("2. Assemble into lipid nanoparticles with NanoAssemblr <https://www.precisionnanosystems.com/platform-technologies/product-comparison/ignite>.")
    print("3. Buffer exchange or something, idk.")
    print("4. Inject into your (least?) favorite mice.")
    print()
    print("For more information, see <https://www.nature.com/articles/s41467-022-30878-4>")
    print()
    print("Thank you for using Autovax 🌈")
    
    
    
if __name__ == "__main__":
    main()



