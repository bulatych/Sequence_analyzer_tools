import os
from typing import Dict, Union, Tuple, List
from additional_modules.dna_rna_tools import (
    seq_transcr,
    seq_reverse,
    seq_compl,
    seq_rev_compl,
    is_valid_sequence,
)
from additional_modules.filter_fastq import (
    is_length_bounds,
    gc_content_calculator,
    calc_quality,
    check_quality,
    check_length,
    check_gc_content,
)


def run_dna_rna_tools(*args: Union[str, List[str]]) -> Union[str, List[str]]:
    """
    Function for working with DNA/RNA sequences

    :param args:  accepts one or more nucleic acid sequences and performs 4
    operations
    :param procedure:
        "transcribe"- transcription of DNA into RNA
        "reverse" - returns the reverse sequence
        "complement" - returns the complementary sequence
        "reverse_complement" - returns the reverse complementary sequence
        The function excludes the input sequence for the simultaneous
        presence of "U", "u"  and "T", "t" in the composition
        Function takes into account the register of nucleotides in the sequence

    The funcion uses additional module dna_rna_tools
    """
    *seqs, procedure = args
    if not is_valid_sequence(seqs):
        return "Invalid sequence"
    operations = {
        "transcribe": seq_transcr,
        "reverse": seq_reverse,
        "complement": seq_compl,
        "reverse_complement": seq_rev_compl,
    }
    if procedure in operations:
        operation = operations[procedure]
    else:
        return "Unknown procedure"

    result = [operation(seq) for seq in seqs]

    return result[0] if len(result) == 1 else result


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: Union[Tuple[float, float], float] = (0, 100),
    length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
    quality_threshold: int = 0,
) -> Dict[str, Tuple[str, str]]:
    """
    The function works with fastq sequences. Accepts 5 arguments as input

    The function works on the fly, accepts a fastq file,
    selects sequences for recording and saves the filtered data. Functions
    checks the existence of output directory and if missing creates
    a folder "filtered"
    :param input_fastq: str
         The path to real fastq file with sequence id, quality, sequence
    :param output_fastq: str
         The name of output fastq file with filtered sequences
    :param gc_bounds:
         The GC interval of the composition for filtering.If there is one
         argument, then it is assumed that this is the
         upper bound
    :param length_bounds:
         Length interval for filtering. It is set similarly to gc_bounds
    :param quality_threshold:
         The threshold value of the average reed quality for filtering.
         Quality is calculated by converting the ASCII encoding scale
         to Q-Score.The average quality is calculated for all nucleotides
         and those below the threshold are discarded.

    Function uses additional module filter_fastq.py
    """
    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, (int)):
        length_bounds = (0, length_bounds)

    output_dir = "filtered"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    output_path = os.path.join(output_dir, output_fastq)

    with open(input_fastq, "r") as file_in, open(output_path, "w") as file_out:
        line_counter = 0
        seq_id, sequence, quality = "", "", ""
        for line in file_in:
            line = line
            if line_counter == 0:
                seq_id = line
            elif line_counter == 1:
                sequence = line
            elif line_counter == 3:
                quality = line

            if (
                check_quality(quality, quality_threshold)
                and check_length(length_bounds, sequence)
                and check_gc_content(gc_bounds, sequence)
            ):
                file_out.write(seq_id)
                file_out.write(sequence)
                file_out.write("+\n")
                file_out.write(quality)
            line_counter = (line_counter + 1) % 4


# Example of using function
filter_fastq("example_fastq.fastq", "filtered_fastq.fastq", (0, 20), (0, 100), 30)
