def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str) -> None:
    """
    Convert a multi-line FASTA file into a single-line FASTA file.

    This function reads a multi-line FASTA formatted file specified by
    `input_fasta`, concatenates the nucleotide or protein sequences into
    single lines, and writes the result to a new FASTA formatted file
    specified by `output_fasta`
    :param input_fasta: str
         The path to the input FASTA file, which can contain
         multiple lines of sequences.
    :param output_fasta: str
         The path to the output FASTA file, where the converted
         single-line sequences will be saved.

    :return: None
         The function does not return any value. It writes to the
         specified output file instead.
    """

    with open(input_fasta, "r") as fasta_in, open(output_fasta, "w") as fasta_out:
        sequence = ""
        for line in fasta_in:
            line = line.strip()
            if line.startswith(">"):
                if sequence:
                    fasta_out.write(sequence + "\n")
                fasta_out.write(line + "\n")
                sequence = ""
            else:
                sequence += line
        if sequence:  # Check if there last sequence remains
            fasta_out.write(sequence + "\n")


def parse_blast_output(input_file: str, output_file: str) -> None:
    """
    Parse a BLAST output file to extract and sort significant protein names.

    This function reads a BLAST output file specified by `input_file`,
    extracts protein names from the section that details sequences producing
    significant alignments, and writes the unique names into a specified
    output file `output_file`. The protein names are sorted in
    alphabetical order before being written.

    :param input_file:
         The path to the input BLAST output file, which contains the results
         of a BLAST search in text format
    :param output_file:
         The path to the output text file where the sorted list of protein
         names will be written.
    :return: None
         This function does not return any value. Instead, it writes the
         sorted protein names directly to the specified output file.
    """
    proteins = []
    with open(input_file, "r") as blast_in:
        for line in blast_in:
            if line.startswith("Sequences producing significant alignments:"):
                blast_in.readline()
                blast_in.readline()
                description_line = blast_in.readline()
                if description_line:
                    description_parts = description_line.split("  ")
                    if description_parts:
                        protein_name = description_parts[0]
                        proteins.append(protein_name)
    sorted_proteins = sorted(proteins, key=str.lower)

    with open(output_file, "w") as file_out:
        for protein in sorted_proteins:
            file_out.write(protein + "\n")
