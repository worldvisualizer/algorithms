# import statements
import time
from collections import defaultdict

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
    dna_gc_data = []
    for dataset in fasta_file.split('>')[1:]:
        split = dataset.split('\n')
        dna_seq = ''.join(split[1:]).strip()
        dna_gc_data.append((split[0], gc_contents(dna_seq)))
    return max(dna_gc_data, key=lambda x: x[1])


if __name__ == '__main__':
    print(gc_contents('AG'))

    