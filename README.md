# Sequence Analyzer

| Sequence Analyzer - A toolkit for processing DNA/RNA sequences and FASTQ files. This tool allows sequence manipulation and filtering based on specific criteria. | <img src="tool_pict.png" alt="Project Logo" width="300"/> |
|:--------------------------------------------------------|:-------------------------------------------------------:|

**Author:**  
**Software, Idea, Testing**: [*Bulat Rakhimov*](https://t.me/bulatych_7)  

## Table of Contents
- [Description](#description)
- [Installation](#installation)
- [Examples](#examples)
- [FAQ](#faq)
- [Contact](#contact)

## Description
The **Sequence Analyzer** toolkit provides the following key features:

- **BiologicalSequence Classes**: Abstract base class and specific implementations for DNA, RNA, and protein sequences, allowing common operations like reverse sequences, transcriptions, complement sequences, and finding stop codons.
  - **DNASequence**: Includes functions for transcription, reverse complement, and sequence validation.
  - **RNASequence**: Implements stop codon identification.
  - **AminoAcidSequence**: Provides functionality for identifying aromatic amino acids and counting residue occurrences.

- **filter_fastq**: Filters FASTQ sequences based on GC content, sequence length, and quality threshold. The function creates an output directory ("filtered") if it does not exist and saves filtered sequences.

## Installation
To use the Sequence Analyzer toolkit, follow these steps:

1. Clone the repository:
   ```bash
   git clone git@github.com:bulatych/Sequence_analyzer_tools.git
   cd Sequence_analyzer_tools
   ```

2. Run the script:
   ```bash
   python Sequence_analyzer.py
   ```

## Examples

### Filtering FASTQ Sequences
The `filter_fastq` function filters sequences based on user-defined criteria:

```python
from Sequence_analyzer import filter_fastq

# Parameters
gc_bounds = (0, 20)
length_bounds = (0, 100)
quality_threshold = 30

filter_fastq("example_fastq.fastq", "filtered_fastq.fastq", gc_bounds, length_bounds, quality_threshold)
```

### DNA/RNA Sequence Operations

```python
from Sequence_analyzer import DNASequence, RNASequence, AminoAcidSequence

# DNA Sequence Example
dna_seq = DNASequence("ATGCGT")
print(dna_seq.reverse())  # Reverse sequence
print(dna_seq.complement())  # Complement sequence
print(dna_seq.transcribe())  # Transcribe to RNA

# RNA Sequence Example
rna_seq = RNASequence("AUGCGU")
print(rna_seq.find_stop_codons(str(rna_seq)))  # Identify stop codons

# Amino Acid Sequence Example
protein_seq = AminoAcidSequence("ACDEFGHIKLM")
print(protein_seq.count_repeating_amino(str(protein_seq)))  # Count amino acids
```

## FAQ

**1. What is the purpose of the Sequence Analyzer toolkit?**  
This toolkit helps researchers process and analyze genomic data efficiently by providing tools for sequence manipulation and filtering FASTQ sequences.

**2. Can I use this toolkit for other sequencing formats?**  
Currently, the toolkit supports DNA, RNA, and FASTQ sequences. Future updates may include additional formats.

**3. Where are the filtered sequences stored?**  
The filtered sequences are saved in the `filtered/` directory created automatically if it does not exist.

## Contact
For any issues or feature requests, feel free to contact: [*Bulat Rakhimov*](https://t.me/bulatych_7)



