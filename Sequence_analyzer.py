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
    quality_check,
)


def run_dna_rna_tools(*args: Union[str, List[str]]) -> Union[str, List[str]]:
    """Function for working with DNA/RNA sequences
    Args:
    args:  accepts one or more nucleic acid sequences and performs 4 operations
    procedure:
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
    seqs: Dict[str, Tuple[str, str]],
    gc_bounds: Union[Tuple[float, float], float] = (0, 100),
    length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
    quality_threshold: int = 0,
) -> Dict[str, Tuple[str, str]]:
    """
    The function works with fastq sequences. Accepts 4 arguments as input:

    Arguments:
    seqs:
    a dictionary consisting of fastq sequences, including the name of
    the reads,sequence and quality
    gc_bounds:
    the GC interval of the composition for filtering.
    If there is one argument, then it is assumed that this is the upper bound
    length_bounds:
    length interval for filtering. It is set similarly to gc_bounds
    quality_threshold:
    the threshold value of the average reed quality for filtering.
    Quality is calculated by converting the ASCII encoding scale to Q-Score.
    The average quality is calculated for all nucleotides and those below the
    threshold are discarded.
    Function uses additional module filter_fastq.py
    Function returns a similar dictionary consisting only of those sequences
    that satisfy all the conditions
    """
    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, (int)):
        length_bounds = (0, length_bounds)

    filtered_seqs = {}

    for name, (sequence, quality) in seqs.items():

        if not is_length_bounds(len(sequence), length_bounds):
            continue
        gc_content = gc_content_calculator(sequence)
        if not gc_bounds[0] <= gc_content <= gc_bounds[1]:
            continue
        if quality_check(quality) < quality_threshold:
            continue

        filtered_seqs[name] = (sequence, quality)
    return filtered_seqs
