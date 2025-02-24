from abc import ABC, abstractmethod


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

    @abstractmethod
    def complement(self):
        pass
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
    def _check_alphabet(self) -> bool:
        amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
        return all(residue in amino_acids for residue in self.sequence)

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
