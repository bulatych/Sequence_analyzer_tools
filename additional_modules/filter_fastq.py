from typing import Tuple, Union, Dict


def gc_content_calculator(sequence):
    length_sequence = len(sequence)
    if length_sequence == 0:
        return 0.0
    gc_count = (
        sequence.count("G")
        + sequence.count("C")
        + sequence.count("g")
        + sequence.count("c")
    )
    if sequence:
        return gc_count / length_sequence * 100
    else:
        return 0


def check_length_bounds(length: int, length_bounds: Tuple[int, int]):
    return length_bounds[0] <= length <= length_bounds[1]


def quality_check(quality):
    if len(quality) == 0:
        return 0.0
    total_quality = 0
    for char in quality:
        total_quality += ord(char) - 33
    return total_quality / len(quality)