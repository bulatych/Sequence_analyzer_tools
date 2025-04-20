from abc import ABC, abstractmethod
import os
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from typing import Dict, Union, Tuple
import argparse
import logging


class BiologicalSequence(ABC):
    """
    An abstract base class representing a biological sequence.
    """
    def __init__(self, sequence: str, alphabet: set):
        """
        Initialize a biological sequence.

        :param sequence: The sequence string.
        :param alphabet: A set of valid characters for the sequence.
        """
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

    @abstractmethod
    def reverse(self):
        pass


class NucleicAcidSequence(BiologicalSequence):
    """
    A base class for nucleic acid sequences (DNA and RNA).
    """
    complement_map = {}  # for polymorphism

    def complement(self):
        if not self.complement_map:
            raise NotImplementedError(
                "complement method must be implemented in subclasses"
            )
        return self.__class__(
            "".join(self.complement_map[base] for base in self.sequence)
        )

    def reverse(self):
        return self.sequence[::-1]

    def reverse_complement(self):
        return self.complement()[::-1]


class DNASequence(NucleicAcidSequence):
    complement_map = {"A": "T", "T": "A", "C": "G", "G": "C"}  # for polymorphism

    def __init__(self, sequence: str):
        super().__init__(sequence, {"A", "T", "C", "G"})

    def transcribe(self):
        return RNASequence(self.sequence.replace("T", "U"))


class RNASequence(NucleicAcidSequence):
    complement_map = {"A": "U", "U": "A", "C": "G", "G": "C"}  # for polymorphism
    stop_codons = {"UAA", "UAG", "UGA"}

    def __init__(self, sequence: str):
        super().__init__(sequence, {"A", "U", "C", "G"})

    def find_stop_codons(self):
        """
        Identify positions of stop codons in an RNA sequence.

        :param rna_sequence: A string representing the RNA sequence.
        :return: A list of indices where stop codons occur.
        """
        stop_positions = []
        for i in range(0, len(self.sequence) - 2, 3):
            codon = self.sequence[i:i + 3]
            if codon in self.stop_codons:
                stop_positions.append(i)
        return stop_positions



class AminoAcidSequence(BiologicalSequence):
    def __init__(self, sequence: str):
        super().__init__(
            sequence,
            {
                "A",
                "C",
                "D",
                "E",
                "F",
                "G",
                "H",
                "I",
                "K",
                "L",
                "M",
                "N",
                "P",
                "Q",
                "R",
                "S",
                "T",
                "V",
                "W",
                "Y",
            },
        )

    def reverse(self):
        return self.sequence[::-1]

    def find_aromatic_acid(self, amino):
        return amino in {"F", "W", "Y"}

    def count_repeating_amino(self, amino_sequence):
        """
        Count occurrences of each amino acid in a sequence.

        :param amino_sequence: A string representing the amino acid sequence.
        :return: A dictionary with amino acids as keys and their counts as values.
        """
        amino_count = {}
        for amino in self.sequence:
            amino_count[amino] = amino_count.get(amino, 0) + 1
        return amino_count



logging.basicConfig(
    filename="filter_fastq.log",  
    level=logging.DEBUG, 
    format="%(asctime)s - %(levelname)s - %(message)s",  
)


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: Union[Tuple[float, float], float] = (0, 100),
    length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
    quality_threshold: int = 0,
) -> Dict[str, Tuple[str, str]]:
    """
    Filter sequences from a FASTQ file based on GC content, length, and quality.

    :param input_fastq: Path to the input FASTQ file.
    :param output_fastq: Path to the output FASTQ file.
    :param gc_bounds: Tuple representing the min and max GC content allowed (percentage).
    :param length_bounds: Tuple representing the min and max length of sequences allowed.
    :param quality_threshold: Minimum average quality score required.
    """
    logging.info(f"Starting to filter the FASTQ file: {input_fastq}")

    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, (int)):
        length_bounds = (0, length_bounds)

    output_path = output_fastq

    try:
        with open(input_fastq, "r") as file_in, open(output_path, "w") as file_out:
            for record in SeqIO.parse(file_in, "fastq"):
                sequence = record.seq
                quality = record.letter_annotations.get("phred_quality", [])
    
                if not (length_bounds[0] <= len(sequence) <= length_bounds[1]):
                    continue
    
                avg_quality = sum(quality) / len(quality) if quality else 0
                if avg_quality < quality_threshold:
                    continue
                gc_content = gc_fraction(sequence) * 100
                if not (gc_bounds[0] <= gc_content <= gc_bounds[1]):
                    continue
    
                SeqIO.write(record, file_out, "fastq")
    
        print(f"All sequences were filtered and saved in: '{output_path}'.")
        
        logging.info(f"Filtering completed. Output file: {output_path}")
        
    except Exception as e:
        logging.error(f"An error occurred: {e}")
        raise
        



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter Fastq sequnces by GC content, quality and length")
    parser.add_argument('-i', "--input", required= True, help = "Path to the input file")
    parser.add_argument('-o', "--output", required= True, help = "Name of the output file")
    parser.add_argument("--gc_bounds", nargs="+", type=float, default=[0, 100], help="GC content bounds.")
    parser.add_argument("--length_bounds", nargs="+", type=int, default=[0, 2**32], help ="Length of bounds")    
    parser.add_argument('-q', "--quality_ths", type=int, default=0, help= "Minimum average quality threshold")
    args = parser.parse_args()
    gc_bounds = tuple(args.gc_bounds) if len(args.gc_bounds) == 2 else args.gc_bounds[0]
    length_bounds = tuple(args.length_bounds) if len(args.length_bounds) == 2 else args.length_bounds[0]

    filter_fastq(
        input_fastq=args.input,
        output_fastq=args.output,
        gc_bounds=gc_bounds,
        length_bounds=length_bounds,
        quality_threshold=args.quality_ths,
    )





