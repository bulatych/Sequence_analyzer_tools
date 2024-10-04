from typing import Dict, Union, Tuple
from additional_modules.dna_rna_tools import seq_transcr, seq_reverse, seq_compl, seq_rev_compl, valid_sequence
def run_dna_rna_tools(*args):
    *seqs, procedure = args
    if not valid_sequence(seqs):
        return "Invalid sequence"
    operations = {
        "transcribe": seq_transcr,
        "reverse": seq_reverse,
        "complement": seq_compl,
        "reverse_complement": seq_rev_compl
    }
    if procedure in operations:
        operation = operations[procedure]
    else:
        return "Unknown procedure"

    result = [operation(seq) for seq in seqs]

    return result[0] if len(result) == 1 else result

from additional_modules.filter_fastq import check_length_bounds, gc_content_calculator, quality_check
def filter_fastq(seqs: Dict[str, Tuple[str, str]],
    gc_bounds: Union[Tuple[float,float], float] = (0, 100),
    length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
    quality_threshold: int = 0) -> Dict[str, Tuple[str, str]]:

    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds,(int)):
        length_bounds = (0, length_bounds)

    filtered_seqs = {}

    for name, (sequence, quality) in seqs.items():

        if not check_length_bounds(len(sequence),length_bounds):
            continue
        gc_content = gc_content_calculator(sequence)
        if not gc_bounds[0] <= gc_content <= gc_bounds[1]:
            continue
        if quality_check(quality) < quality_threshold:
            continue

        filtered_seqs[name] = (sequence, quality)
    return filtered_seqs
