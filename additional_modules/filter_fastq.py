from tabnanny import check
from typing import Tuple


def gc_content_calculator(sequence: str) -> float:
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


def is_length_bounds(length: int, length_bounds: Tuple[int, int]) -> bool:
    return length_bounds[0] <= length <= length_bounds[1]


def calc_quality(quality: str) -> float:
    if len(quality) == 0:
        return 0.0
    total_quality = 0
    for char in quality:
        total_quality += ord(char) - 33
    return total_quality / len(quality)

def quality_check(quality: str, quality_threshold: int) -> bool:
    avg_quality = calc_quality(quality)
    return avg_quality >= quality_threshold

def check_length(length_bounds: Tuple[int, int], sequence: str) -> bool:
    return is_length_bounds(len(sequence), length_bounds)


def check_gc_content(gc_bounds: Tuple[float, float], sequence: str) -> bool:
    gc_content = gc_content_calculator(sequence)
    return gc_bounds[0] <= gc_content <= gc_bounds[1]