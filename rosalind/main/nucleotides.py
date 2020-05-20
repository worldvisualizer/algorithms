# import statements
import time
from collections import defaultdict

CODON_MAP = {
    "UUU": "F",
    "CUU": "L",
    "AUU": "I",
    "GUU": "V",
    "UUC": "F",
    "CUC": "L",
    "AUC": "I",
    "GUC": "V",
    "UUA": "L",
    "CUA": "L",
    "AUA": "I",
    "GUA": "V",
    "UUG": "L",
    "CUG": "L",
    "AUG": "M",
    "GUG": "V",
    "UCU": "S",
    "CCU": "P",
    "ACU": "T",
    "GCU": "A",
    "UCC": "S",
    "CCC": "P",
    "ACC": "T",
    "GCC": "A",
    "UCA": "S",
    "CCA": "P",
    "ACA": "T",
    "GCA": "A",
    "UCG": "S",
    "CCG": "P",
    "ACG": "T",
    "GCG": "A",
    "UAU": "Y",
    "CAU": "H",
    "AAU": "N",
    "GAU": "D",
    "UAC": "Y",
    "CAC": "H",
    "AAC": "N",
    "GAC": "D",
    "UAA": "Stop",
    "CAA": "Q",
    "AAA": "K",
    "GAA": "E",
    "UAG": "Stop",
    "CAG": "Q",
    "AAG": "K",
    "GAG": "E",
    "UGU": "C",
    "CGU": "R",
    "AGU": "S",
    "GGU": "G",
    "UGC": "C",
    "CGC": "R",
    "AGC": "S",
    "GGC": "G",
    "UGA": "Stop",
    "CGA": "R",
    "AGA": "R",
    "GGA": "G",
    "UGG": "W",
    "CGG": "R",
    "AGG": "R",
    "GGG": "G"
}


def count_atcgs(dna_sequence):
    """
        Parameters
        ----------
        dna_sequence : str
            Sequence of ATCG

        Returns
        -------
        defaultdict(int)
            a dict of nucleotide and count
    """
    nuc_counts = defaultdict(int)
    for nuc in dna_sequence:
        nuc_counts[nuc] += 1
    return nuc_counts


def dna_transcription(dna_sequence):
    # reference: https://waymoot.org/home/python_string/
    """
        Parameters
        ----------
        dna_sequence : str
            Sequence of ATCG

        Returns
        -------
        str
            Sequence of AUCG
    """
    return ''.join([seq if seq != 'T' else 'U' for seq in dna_sequence])


def rna_translation(rna_sequence):
    """
        Parameters
        ----------
        rna_sequence : str
            Sequence of AUCG

        Returns
        -------
        str
            Amino acid sequence
    """
    amino_acids = []
    rna_array = list(rna_sequence)
    codons = [
        ''.join(rna_array[i:i+3])
        for i in range(0, len(rna_array), 3)
    ]
    for cod in codons:
        aa = CODON_MAP[cod]
        if aa == 'Stop':
            break
        amino_acids.append(aa)
    return ''.join(amino_acids)


def reverse_complement(dna_sequence):
    # reference: https://stackoverflow.com/questions/36949665/fastest-way-to-reverse-a-string-in-python
    # reference: https://stackoverflow.com/questions/1228299/changing-one-character-in-a-string-in-python/22149018#22149018
    """
        Parameters
        ----------
        dna_sequence : str
            Sequence of ATCG, must be ascii compliant

        Returns
        -------
        str
            Reverse + Complement of dna_sequence, ascii
    """
    complement_map = {
        65: 84,
        84: 65,
        67: 71,
        71: 67
    }
    dna_array = bytearray(dna_sequence, 'ascii')
    for index in range(len(dna_array)):
        dna_array[index] = complement_map[dna_array[index]]
    return dna_array[::-1].decode('ascii')


def gc_contents(dna_sequence):
    # Python's float type is C double in CPython
    """
        Parameters
        ----------
        dna_sequence : str
            Sequence of ATCG, must be ascii compliant

        Returns
        -------
        float
            GC Contents of DNA sequence
    """
    count = 0
    for nuc in dna_sequence:
        if nuc in ['G', 'C']:
            count += 1
    return float(count) / len(dna_sequence) * 100


def top_gc_contents_fasta(fasta_file):
    """
        Parameters
        ----------
        fasta_file : str
            fasta contents in a single string.
            > denotes metadata
            rest of the lines are DNA contents.
            e.g.)

        >Rosalind_6404
        CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCC
        TCCCACTAATAATTCTGAGG
        >Rosalind_5959
        CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCT
        ATATCCATTTGTCAGCAGACACGC

        Returns
        -------
        tuple(metadata, gc_contents)
            metadata : str
                metadata about the DNA
            gc_contents : float
                GC Contents of the DNA sequence
    """
    dna_gc_data = []
    for dataset in fasta_file.split('>')[1:]:
        split = dataset.split('\n')
        dna_seq = ''.join(split[1:]).strip()
        dna_gc_data.append((split[0], gc_contents(dna_seq)))
    return max(dna_gc_data, key=lambda x: x[1])


if __name__ == '__main__':
    pass    