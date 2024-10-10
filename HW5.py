import os
from typing import Dict, Union, Tuple, List
from additional_modules.filter_fastq import (
    is_length_bounds,
    gc_content_calculator,
    calc_quality,
    quality_check,
    check_length,
    check_gc_content
)
def filter_fastq(
    input_fastq, output_fastq,
    gc_bounds: Union[Tuple[float, float], float] = (0, 100),
    length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
    quality_threshold: int = 0,) -> Dict[str, Tuple[str, str]]:
    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, (int)):
        length_bounds = (0, length_bounds)


    output_dir = "filtered"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    output_path = os.path.join(output_dir, output_fastq)

    with (open(input_fastq, "r") as file_in, open(output_path, "w") as file_out):
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

            if (quality_check(quality, quality_threshold) and
               check_length(length_bounds, sequence) and
               check_gc_content(gc_bounds, sequence)):
                    file_out.write(seq_id)
                    file_out.write(sequence)
                    file_out.write("+\n")
                    file_out.write(quality)
            line_counter = (line_counter + 1) % 4
filter_fastq("example_fastq.fastq", "filtered_fastq.fastq", (0, 20), (0, 100), 30)
