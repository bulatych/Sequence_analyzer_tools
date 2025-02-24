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
        """Проверка, что все символы последовательности присутствуют в алфавите."""
        return all(base in self.alphabet for base in self.sequence)

    @abstractmethod
    def complement(self):
        pass


class NucleicAcidSequence(BiologicalSequence):
    complement_map = {}

    def complement(self):
        if not self.complement_map:
            raise NotImplementedError("complement method must be implemented in subclasses")
        return self.__class__("".join(self.complement_map[base] for base in self.sequence))

    def reverse(self):
        return self.sequence[::-1]

    def reverse_complement(self):
        return self.complement()[::-1]


class DNASequence(NucleicAcidSequence):
    complement_map = {"A": "T", "T": "A", "C": "G", "G": "C"}

    def __init__(self, sequence: str):
        super().__init__(sequence, {"A", "T", "C", "G"})

    def transcribe(self):
        return RNASequence(self.sequence.replace("T", "U"))


class RNASequence(NucleicAcidSequence):
    complement_map = {"A": "U", "U": "A", "C": "G", "G": "C"}

    def __init__(self, sequence: str):
        super().__init__(sequence, {"A", "U", "C", "G"})


class AminoAcidSequence(BiologicalSequence):
    def __init__(self, sequence: str):
        # Предположим, что алфавит аминокислот - это стандартный набор 20 аминокислот
        amino_acid_alphabet = set("ACDEFGHIKLMNPQRSTVWY")
        super().__init__(sequence, amino_acid_alphabet)


# class AminoAcidSequence(BiologicalSequence):
#     def __init__(self,sequence):
#         super().__init__(sequence)

rna = RNASequence('AGCA')
rna.reverse_complement()
print(rna.reverse_complement())
dna = DNASequence('AT')
print(dna.reverse())
