# Codon table (There must be an easier way rather than writing all of this manually...)
codon_table = {
    'UUU': 'Phe (F)', 'UUC': 'Phe (F)', 'UUA': 'Leu (L)', 'UUG': 'Leu (L)',
    'CUU': 'Leu (L)', 'CUC': 'Leu (L)', 'CUA': 'Leu (L)', 'CUG': 'Leu (L)',
    'AUU': 'Ile (I)', 'AUC': 'Ile (I)', 'AUA': 'Ile (I)', 'AUG': 'Met (M)',
    'GUU': 'Val (V)', 'GUC': 'Val (V)', 'GUA': 'Val (V)', 'GUG': 'Val (V)',
    'UCU': 'Ser (S)', 'UCC': 'Ser (S)', 'UCA': 'Ser (S)', 'UCG': 'Ser (S)',
    'CCU': 'Pro (P)', 'CCC': 'Pro (P)', 'CCA': 'Pro (P)', 'CCG': 'Pro (P)',
    'ACU': 'Thr (T)', 'ACC': 'Thr (T)', 'ACA': 'Thr (T)', 'ACG': 'Thr (T)',
    'GCU': 'Ala (A)', 'GCC': 'Ala (A)', 'GCA': 'Ala (A)', 'GCG': 'Ala (A)',
    'UAU': 'Tyr (Y)', 'UAC': 'Tyr (Y)', 'UAA': 'STOP', 'UAG': 'STOP',
    'CAU': 'His (H)', 'CAC': 'His (H)', 'CAA': 'Gln (Q)', 'CAG': 'Gln (Q)',
    'AAU': 'Asn (N)', 'AAC': 'Asn (N)', 'AAA': 'Lys (K)', 'AAG': 'Lys (K)',
    'GAU': 'Asp (D)', 'GAC': 'Asp (D)', 'GAA': 'Glu (E)', 'GAG': 'Glu (E)',
    'UGU': 'Cys (C)', 'UGC': 'Cys (C)', 'UGA': 'STOP', 'UGG': 'Trp (W)',
    'CGU': 'Arg (R)', 'CGC': 'Arg (R)', 'CGA': 'Arg (R)', 'CGG': 'Arg (R)',
    'AGU': 'Ser (S)', 'AGC': 'Ser (S)', 'AGA': 'Arg (R)', 'AGG': 'Arg (R)',
    'GGU': 'Gly (G)', 'GGC': 'Gly (G)', 'GGA': 'Gly (G)', 'GGG': 'Gly (G)'
}


def complement_dna(dna_seq):
    """Generate the complement of the DNA strand"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join([complement[base] for base in dna_seq])


def transcribe_to_mrna(dna_seq):
    """Transcribe DNA to mRNA"""
    return dna_seq.replace('T', 'U')


def translate_to_amino_acids(mrna_seq):
    """Translate mRNA to amino acid seq"""
    amino_acids = []
    for i in range(0, len(mrna_seq), 3):
        codon = mrna_seq[i:i + 3]
        amino_acid = codon_table.get(codon, 'Unknown')
        if amino_acid == 'STOP':  # Translation will stop at the STOP codon
            break
        amino_acids.append(amino_acid)
    return ' – '.join(amino_acids)


def translate_dna_to_protein(dna_seq):
    if len(dna_seq) % 3 != 0:
        raise ValueError("The input DNA seq length must be a multiple of 3.")

    print(f"Input DNA = {dna_seq}")

    complement = complement_dna(dna_seq)
    print(f"Complement = {complement}")

    mrna = transcribe_to_mrna(complement)
    print(f"mRNA = {mrna}")

    amino_acid_seq = translate_to_amino_acids(mrna)
    print(f"Aminoacid = {amino_acid_seq}")


# Example usage:
translate_dna_to_protein("TTACGA")
translate_dna_to_protein("AAACGA")
