# Dictionary mapping amino acids to codons (Same thing as before, must be an easier way than doing this manually...)
codon_table = {
    'F': ['UUU', 'UUC'],                                # Phenylalanine
    'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],    # Leucine
    'I': ['AUU', 'AUC', 'AUA'],                         # Isoleucine
    'M': ['AUG'],                                       # Methionine (Start Codon)
    'V': ['GUU', 'GUC', 'GUA', 'GUG'],                  # Valine
    'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],    # Serine
    'P': ['CCU', 'CCC', 'CCA', 'CCG'],                  # Proline
    'T': ['ACU', 'ACC', 'ACA', 'ACG'],                  # Threonine
    'A': ['GCU', 'GCC', 'GCA', 'GCG'],                  # Alanine
    'Y': ['UAU', 'UAC'],                                # Tyrosine
    'H': ['CAU', 'CAC'],                                # Histidine
    'Q': ['CAA', 'CAG'],                                # Glutamine
    'N': ['AAU', 'AAC'],                                # Asparagine
    'K': ['AAA', 'AAG'],                                # Lysine
    'D': ['GAU', 'GAC'],                                # Aspartic acid
    'E': ['GAA', 'GAG'],                                # Glutamic acid
    'C': ['UGU', 'UGC'],                                # Cysteine
    'W': ['UGG'],                                       # Tryptophan
    'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],    # Arginine
    'G': ['GGU', 'GGC', 'GGA', 'GGG'],                  # Glycine
    '*': ['UAA', 'UAG', 'UGA'],                         # Stop codons
}


def calculate_codon_frequency(mrna_seq, amino_acid_seq):
    codon_count = {}

    for amino_acid in amino_acid_seq:
        if amino_acid in codon_table:
            possible_codons = codon_table[amino_acid]
            for codon in possible_codons:
                codon_count[codon] = 0

            for i in range(0, len(mrna_seq) - 2, 3):
                codon = mrna_seq[i:i + 3]
                if codon in possible_codons:
                    codon_count[codon] += 1

    return codon_count


input_amino_acid_seq = input("Input amino acid seq (max 3): ").upper()
if len(input_amino_acid_seq) > 3:
    print("Please enter a maximum of 3 amino acids.")

input_mrna_seq = input("Input mRNA seq: ").upper()

codon_frequencies = calculate_codon_frequency(input_mrna_seq, input_amino_acid_seq)

print(f"mRNA = {input_mrna_seq}")
for result_codon, count in codon_frequencies.items():
    if count != 0:
        print(f"{result_codon} = {count}")

