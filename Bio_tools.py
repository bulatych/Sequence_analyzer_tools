from abc import ABC, abstractmethod
import os
from os import write
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction
from typing import Dict, Union, Tuple



class BiologicalSequence(ABC):
    def __init__(self, sequence: str, alphabet: set):
        self.sequence = sequence
        self.alphabet = alphabet
        if not self._check_alphabet():
            raise ValueError("Invalid sequence alphabet")

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, index):
        return self.sequence[index]

    def __str__(self):
        return self.sequence

    def __repr__(self):
        return f"{self.__class__.__name__}('{self.sequence}')"

    def _check_alphabet(self) -> bool:
        return all(base in self.alphabet for base in self.sequence)

    # @abstractmethod
    # def complement(self):
    #     pass
    @abstractmethod
    def reverse(self):
        pass

class NucleicAcidSequence(BiologicalSequence):
    complement_map = {} # for polymorphism

    def complement(self):
        if not self.complement_map:
            raise NotImplementedError("complement method must be implemented in subclasses")
        return self.__class__("".join(self.complement_map[base] for base in self.sequence))

    def reverse(self):
        return self.sequence[::-1]

    def reverse_complement(self):
        return self.complement()[::-1]


class DNASequence(NucleicAcidSequence):
    complement_map = {"A": "T", "T": "A", "C": "G", "G": "C"} # for polymorphism

    def __init__(self, sequence: str):
        super().__init__(sequence, {"A", "T", "C", "G"})

    def transcribe(self):
        return RNASequence(self.sequence.replace("T", "U"))


class RNASequence(NucleicAcidSequence):
    complement_map = {"A": "U", "U": "A", "C": "G", "G": "C"} # for polymorphism
    stop_codons = {"UAA", "UAG", "UGA"}
    def __init__(self, sequence: str):
        super().__init__(sequence, {"A", "U", "C", "G"})

    def find_stop_codons(self, rna_sequence: str):
        stop_positions = []
        for i in range(0, len(rna_sequence) - 2, 3):
            codon = rna_sequence[i:i + 3]
            if codon in self.stop_codons:
                stop_positions.append(i)
        return stop_positions


class AminoAcidSequence(BiologicalSequence):
    def __init__(self, sequence: str):
        super().__init__(sequence, { 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N',
                                   'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'
})
    def reverse(self):
        return self.sequence[::-1]

    def find_aromatic_acid(self, amino):
        return amino in {'F', 'W', 'Y'}

    def count_repeating_amino(self, amino_sequence):
        amino_count = {}
        for amino in amino_sequence:
            if amino in amino_count:
                amino_count[amino] += 1
            else:
                amino_count[amino] = 1
        return amino_count

    def _check_alphabet(self) -> bool:
        amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
        return all(residue in amino_acids for residue in self.sequence)



def filter_fastq(input_fastq: str,
    output_fastq: str,
    gc_bounds: Union[Tuple[float, float], float] = (0, 100),
    length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
    quality_threshold: int = 0,
) -> Dict[str, Tuple[str, str]]:

    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, (int)):
        length_bounds = (0, length_bounds)

    output_dir = "filtered"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    output_path = os.path.join(output_dir, output_fastq)

    with open(input_fastq, 'r') as file_in, open(output_path, 'w') as file_out:
        for record in SeqIO.parse(file_in, 'fastq'):
            sequence = record.seq
            quality = record.letter_annotations.get('phred_quality', [])

            if not (length_bounds[0] <= len(sequence) <= length_bounds[1]):
                continue

            avg_quality = sum(quality) / len(quality) if quality else 0
            if avg_quality < quality_threshold:
                continue
            gc_content = gc_fraction(sequence) * 100
            if not (gc_bounds[0] <= gc_content <= gc_bounds[1]):
                continue

            SeqIO.write(record, file_out, 'fastq')

    print(f"All sequences were filtered and saved in: '{output_path}'.")

