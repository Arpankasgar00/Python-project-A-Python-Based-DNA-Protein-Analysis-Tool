"""
=============================================================
  BioSeq Toolkit: A Python-Based DNA & Protein Analysis Tool
=============================================================
  Author  : Arpan kasgar
  Purpose : Analyze DNA sequences (input, validate, analyze,
    translate, find motifs, detect mutations, FASTA)

  
"""

# ─────────────────────────────────────────────
#  CODON TABLE  (RNA codon → Amino Acid)
# ─────────────────────────────────────────────
CODON_TABLE = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


# ═══════════════════════════════════════════════════════════
#  STEP 1 — INPUT SYSTEM
# ═══════════════════════════════════════════════════════════

def get_dna_input():
    """
    Step 1: Accept DNA from user (manual typing OR FASTA file).
    Returns the cleaned, uppercase DNA string.

    """
    print("\n" + "="*55)
    print("  STEP 1 — DNA INPUT")
    print("="*55)
    print("  [1] Type / paste sequence manually")
    print("  [2] Load from a FASTA file")
    print("-"*55)

   
    while True:
        choice = input("  Your choice (1 or 2): ").strip()

        if choice == "1":
            dna = input("  Enter DNA sequence: ").strip()
            break

        elif choice == "2":
            filename = input("  Enter FASTA file path: ").strip()
            dna = ""
            try:
                with open(filename, "r") as file:
                    for line in file:
                        if line.startswith(">"):  
                            continue
                        dna += line.strip()       
                break
            except FileNotFoundError:
                print(f"\n  File '{filename}' not found.")
                print("  Please check the path and try again.\n")
                # Loop continues — user re-enters choice

        else:
            print("  Invalid choice. Please enter 1 or 2.\n")
            # Loop continues — user re-enters choice

    # ── Clean the sequence ──
    dna = dna.upper().replace(" ", "")
    print(f"\n  Processed sequence : {dna}")
    return dna


# ═══════════════════════════════════════════════════════════
#  STEP 2 — VALIDATION
# ═══════════════════════════════════════════════════════════

def validate_dna(dna):
    """
    Step 2: Check that every character is A, T, G, or C.
    Returns True if valid, False otherwise.
    """
    print("\n" + "="*55)
    print("  STEP 2 — DNA VALIDATION")
    print("="*55)

    if not dna:
        print("  Sequence is empty — nothing to validate.")
        return False

    valid_bases = "ATGC"
    is_valid = True

    for base in dna:
        if base not in valid_bases:
            is_valid = False
            print(f"  Invalid character found: '{base}'")
            break

    if is_valid:
        print("  Sequence is valid DNA (only A, T, G, C found).")
    else:
        print("  Sequence is invalid — please check your input.")

    return is_valid




# ═══════════════════════════════════════════════════════════
#  STEP 3 — DNA ANALYSIS  (base counts + GC content)
# ═══════════════════════════════════════════════════════════

def analyze_dna(dna):
    """
    Step 3: Count each base, calculate length and GC content.
    Returns a dict with all statistics.
    """
    print("\n" + "="*55)
    print("  STEP 3 — DNA ANALYSIS")
    print("="*55)

    count_a = count_t = count_g = count_c = 0

    for base in dna:
        if base == "A":
            count_a += 1
        elif base == "T":
            count_t += 1
        elif base == "G":
            count_g += 1
        elif base == "C":
            count_c += 1

    length = len(dna)
    gc_content = ((count_g + count_c) / length * 100) if length > 0 else 0.0

    print(f"  Length    : {length}")
    print(f"  A         : {count_a}")
    print(f"  T         : {count_t}")
    print(f"  G         : {count_g}")
    print(f"  C         : {count_c}")
    print(f"  GC Content: {gc_content:.2f}%")

    return {
        "length": length,
        "A": count_a, "T": count_t, "G": count_g, "C": count_c,
        "gc_content": round(gc_content, 2),
    }




# ═══════════════════════════════════════════════════════════
#  STEP 4 — REVERSE COMPLEMENT
# ═══════════════════════════════════════════════════════════

def reverse_complement(dna):
    """
    Step 4: Compute the reverse complement of the DNA strand.
    Complement: A <-> T, G <-> C, then reverse.
    Returns the reverse-complement string.
    """
    print("\n" + "="*55)
    print("  STEP 4 — REVERSE COMPLEMENT")
    print("="*55)

    complement_map = {"A": "T", "T": "A", "G": "C", "C": "G"}

    complement = ""
    for base in dna:
        complement += complement_map[base]

    rev_comp = complement[::-1]   # reverse using slicing

    print(f"  Original          : {dna}")
    print(f"  Complement        : {complement}")
    print(f"  Reverse Complement: {rev_comp}")

    return rev_comp



# ═══════════════════════════════════════════════════════════
#  STEP 5 — TRANSCRIPTION + TRANSLATION (DNA → RNA → Protein)
# ═══════════════════════════════════════════════════════════

def translate_dna(dna):
    """
    Step 5:
      1. Transcribe DNA -> RNA  (T -> U)
      2. Read RNA in codons of 3
      3. Convert codons -> amino acids using CODON_TABLE
      4. Stop at STOP codon (*)
    Returns (rna_string, protein_string).

    """
    print("\n" + "="*55)
    print("  STEP 5 — TRANSCRIPTION & TRANSLATION")
    print("="*55)

    # ── Transcription: DNA -> RNA ──
    rna = dna.replace("T", "U")
    print(f"  DNA : {dna}")
    print(f"  RNA : {rna}")

    # ── Translation: RNA -> Protein ──
    protein = ""
    codon_list = []

   
    for i in range(0, len(rna), 3):
        codon = rna[i:i+3]
        if len(codon) < 3:           # incomplete codon at end — skip
            break
        codon_list.append(codon)
        amino_acid = CODON_TABLE.get(codon, "?")
        if amino_acid == "*":        # STOP codon
            codon_list[-1] = codon + "(STOP)"
            break
        protein += amino_acid

    print(f"  Codons  : {' | '.join(codon_list)}")
    print(f"  Protein : {protein if protein else '(no valid protein found)'}")
    return rna, protein




# ═══════════════════════════════════════════════════════════
#  STEP 6 — MOTIF FINDER  (pattern search)
# ═══════════════════════════════════════════════════════════

def find_motif(dna):
    """
    Step 6: Ask user for a motif pattern, then find all
    positions where it occurs in the DNA sequence.
    Returns (motif, positions_list).

   
       """
    print("\n" + "="*55)
    print("  STEP 6 — MOTIF FINDER")
    print("="*55)

    motif = input("  Enter motif to search (e.g., ATG): ").strip().upper()

    
    if not motif:
        print("  No motif entered — skipping motif search.")
        return "", []

    motif_len = len(motif)
    positions = []

    for i in range(len(dna) - motif_len + 1):
        if dna[i:i + motif_len] == motif:
            positions.append(i)

    if positions:
        print(f"  Motif '{motif}' found {len(positions)} time(s).")
        print(f"  Positions (0-indexed): {positions}")
    else:
        print(f"  Motif '{motif}' not found in the sequence.")

    return motif, positions



# ═══════════════════════════════════════════════════════════
#  STEP 7 — MUTATION DETECTION
# ═══════════════════════════════════════════════════════════

def detect_mutations(dna1):
    """
    Step 7: Accept a second DNA sequence, compare it
    position-by-position with dna1, and report mismatches.
    Returns list of (position, original_base, mutated_base) tuples.

   
    """
    print("\n" + "="*55)
    print("  STEP 7 — MUTATION DETECTION")
    print("="*55)

    dna2 = input("  Enter second DNA sequence to compare: ").strip().upper().replace(" ", "")

   
    if not dna2:
        print("  No sequence entered — skipping mutation detection.")
        return []

    mutations = []

    if len(dna1) != len(dna2):
        print(f"  Sequences differ in length ({len(dna1)} vs {len(dna2)}).")
        print("  Comparing up to the length of the shorter sequence ...")

    compare_len = min(len(dna1), len(dna2))

    for i in range(compare_len):
        if dna1[i] != dna2[i]:
            mutations.append((i, dna1[i], dna2[i]))

    if mutations:
        print(f"\n  Total mutations found: {len(mutations)}")
        for pos, orig, mut in mutations:
            print(f"    Position {pos}: {orig} -> {mut}")
    else:
        print("  No mutations found — sequences are identical.")

    return mutations


# ═══════════════════════════════════════════════════════════
#  STEP 8 — FASTA FILE: MULTIPLE SEQUENCES
# ═══════════════════════════════════════════════════════════

def process_fasta_file():
    """
    Step 8: Read a FASTA file with multiple sequences,
    extract each sequence ID + sequence, then compute
    length and GC content for every entry.
    Returns a list of dicts: [{id, sequence, length, gc_content}, ...]

   
    """
    print("\n" + "="*55)
    print("  STEP 8 — FASTA MULTI-SEQUENCE ANALYSIS")
    print("="*55)

    filename = input("  Enter FASTA file path: ").strip()
    sequences = []
    current_id = None
    current_seq = ""

    try:
        with open(filename, "r") as f:
            for line in f:
                line = line.strip()

                if not line:          
                    continue

                if line.startswith(">"):
                    # Save the previous sequence (if any)
                    if current_id is not None:
                        sequences.append({"id": current_id, "sequence": current_seq})
                    current_id = line[1:]   # strip the '>'
                    current_seq = ""
                else:
                    current_seq += line.upper()

            # Don't forget the last sequence
            if current_id is not None:
                sequences.append({"id": current_id, "sequence": current_seq})

    except FileNotFoundError:
        print(f"  File '{filename}' not found.")
        return []

    if not sequences:
        print("  No sequences found in file.")
        return []

    print(f"\n  Found {len(sequences)} sequence(s):\n")

    for seq_data in sequences:
        seq = seq_data["sequence"]
        length = len(seq)
        gc = ((seq.count("G") + seq.count("C")) / length * 100) if length > 0 else 0.0
        seq_data["length"] = length
        seq_data["gc_content"] = round(gc, 2)

        print(f"  Sequence : {seq_data['id']}")
        print(f"  Length   : {length}")
        print(f"  GC Content: {gc:.2f}%")
        print("  " + "-"*40)

    return sequences


# ═══════════════════════════════════════════════════════════
#  SAVE RESULTS TO FILE
# ═══════════════════════════════════════════════════════════

def save_results(dna, analysis, rev_comp, rna, protein, motif, positions, mutations):
  
  
    """
    Save all analysis results to result.txt

    """


    with open("result.txt", "w", encoding="utf-8") as f:   
        f.write("=" * 55 + "\n")
        f.write("  BioSeq Toolkit v1.1 -- Analysis Report\n")
        f.write("=" * 55 + "\n\n")

        f.write("-- INPUT SEQUENCE --\n")
        f.write(f"DNA: {dna}\n\n")

        f.write("-- DNA ANALYSIS --\n")
        f.write(f"Length    : {analysis['length']}\n")
        f.write(f"A: {analysis['A']}  T: {analysis['T']}  "
                f"G: {analysis['G']}  C: {analysis['C']}\n")
        f.write(f"GC Content: {analysis['gc_content']}%\n\n")

        f.write("-- REVERSE COMPLEMENT --\n")
        f.write(f"{rev_comp}\n\n")

        f.write("-- TRANSCRIPTION & TRANSLATION --\n")
        f.write(f"RNA    : {rna}\n")
        f.write(f"Protein: {protein if protein else '(none)'}\n\n")

        f.write("-- MOTIF SEARCH --\n")
        if motif:
            if positions:
                f.write(f"Motif '{motif}' found {len(positions)} time(s) at: {positions}\n\n")
            else:
                f.write(f"Motif '{motif}' not found.\n\n")
        else:
            f.write("No motif searched.\n\n")

        f.write("-- MUTATION DETECTION --\n")
        if mutations:
            f.write(f"Total mutations: {len(mutations)}\n")
            for pos, orig, mut in mutations:
                f.write(f"  Position {pos}: {orig} -> {mut}\n")
        else:
            f.write("No mutations detected.\n")

    print("\n  Results saved to result.txt")


# ═══════════════════════════════════════════════════════════
#  MAIN DRIVER
# ═══════════════════════════════════════════════════════════

def main():
    print("\n" + "="*55)
    print("  BioSeq Toolkit v1.1")
    print("  DNA & Protein Analysis Tool")
    print("="*55)

    # Step 1: Input
    dna = get_dna_input()

    # Step 2: Validate
    is_valid = validate_dna(dna)
    if not is_valid:
        print("\n  Exiting — fix the sequence and try again.\n")
        return

    # Step 3: Analyze
    analysis = analyze_dna(dna)

    # Step 4: Reverse complement
    rev_comp = reverse_complement(dna)

    # Step 5: Translate
    rna, protein = translate_dna(dna)

    # Step 6: Motif finder
    motif, positions = find_motif(dna)

    # Step 7: Mutation detection
    mutations = detect_mutations(dna)

    # Save output
    save_results(dna, analysis, rev_comp, rna, protein, motif, positions, mutations)

    print("\n" + "="*55)
    print("  Analysis complete!")
    print("  Run Step 8 separately if you need multi-FASTA analysis.")
    print("="*55 + "\n")

    # Step 8: Optional multi-FASTA
    run_fasta = input("  Run multi-FASTA analysis now? (y/n): ").strip().lower()
    if run_fasta == "y":
        process_fasta_file()


if __name__ == "__main__":
    main()
