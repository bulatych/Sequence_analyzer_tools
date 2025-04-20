import unittest
import os
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio_tools import filter_fastq
from unittest.mock import patch

class TestFilterFastq(unittest.TestCase):

    def test_empty_input_file(self): # test for an empty input file
        """
        Checking that the filtering works with an empty input file.
        """
        input_fastq = "test_data/empty.fastq"
        output_fastq = "test_data/output_empty.fastq"
        with open(input_fastq, 'w') as file:
            pass
        filter_fastq(input_fastq, output_fastq, (40, 60), (100, 500), 30)
        self.assertTrue(os.path.exists(output_fastq))
        with open(output_fastq, 'r') as file:
            records = list(SeqIO.parse(file, "fastq"))
            self.assertEqual(len(records), 0)
        os.remove(input_fastq)
        os.remove(output_fastq)

    def test_gc_filter(self): #  test the filtering by GC content
        """
        Checking filtering by GC content.
        """
        input_fastq = "test_data/input_gc.fastq"
        output_fastq = "test_data/output_gc.fastq"
        filter_fastq(input_fastq, output_fastq, (50, 60), (100, 500), 30)
        self.assertTrue(os.path.exists(output_fastq))
        with open(output_fastq, 'r') as file:
            for record in SeqIO.parse(file, "fastq"):
                gc = gc_fraction(record.seq) * 100
                self.assertGreaterEqual(gc, 50)
                self.assertLessEqual(gc, 60)
        os.remove(output_fastq)

    def test_length_filter(self): # test the filtering by sequence length
        """
        Checking filtering by sequence length
        """
        input_fastq = "test_data/input_length.fastq"
        output_fastq = "test_data/output_length.fastq"
        filter_fastq(input_fastq, output_fastq, (0, 100), (150, 200), 30)
        self.assertTrue(os.path.exists(output_fastq))
        with open(output_fastq, 'r') as file:
            for record in SeqIO.parse(file, "fastq"):
                self.assertGreaterEqual(len(record), 150)
                self.assertLessEqual(len(record), 200)
        os.remove(output_fastq)

    def test_quality_filter(self): # test the filtering by sequence quality
        """
        Checking filtering by sequence quality
        """
        input_fastq = "test_data/input_quality.fastq"
        output_fastq = "test_data/output_quality.fastq"
        filter_fastq(input_fastq, output_fastq, (0, 100), (0, 500), 35)
        self.assertTrue(os.path.exists(output_fastq))
        with open(output_fastq, 'r') as file:
            for record in SeqIO.parse(file, "fastq"):
                qual = record.letter_annotations["phred_quality"]
                avg = sum(qual) / len(qual)
                self.assertGreaterEqual(avg, 35)
        os.remove(output_fastq)

    def test_combined_filter(self): # test the combination of all filters
        """
        Checking combination of all filters
        """
        input_fastq = "test_data/input_combined.fastq"
        output_fastq = "test_data/output_combined.fastq"
        filter_fastq(input_fastq, output_fastq, (40, 60), (150, 300), 30)
        self.assertTrue(os.path.exists(output_fastq))
        with open(output_fastq, 'r') as file:
            for record in SeqIO.parse(file, "fastq"):
                gc = gc_fraction(record.seq) * 100
                self.assertGreaterEqual(gc, 40)
                self.assertLessEqual(gc, 60)
                self.assertGreaterEqual(len(record), 150)
                self.assertLessEqual(len(record), 300)
                qual = record.letter_annotations["phred_quality"]
                avg = sum(qual) / len(qual)
                self.assertGreaterEqual(avg, 30)
        os.remove(output_fastq)

    @patch("logging.error")
    def test_logging_error_invalid_input(self, mock_error): # test logging check for incorrect parameters
        """
        Checking that an error is logged with incorrect input data.
        """
        input_fastq = "test_data/invalid_input.fastq"
        output_fastq = "test_data/output_error.fastq"
        filter_fastq(input_fastq, output_fastq, (40, 60), (100, 500), 30)
        mock_error.assert_called_with("Файл не найден: test_data/invalid_input.fastq")
        
    def test_output_file_created(tmp_path):
        """
        Checking output file is being created.
        """
        input_fastq = "test_data/sample.fastq"
        output_fastq = "filtered_result.fastq"
        output_path = tmp_path / "filtered" / output_fastq
        filter_fastq(input_fastq, output_fastq, gc_bounds=(0, 100), length_bounds=(0, 1000), quality_threshold=0)
        assert output_path.exists()


    @patch("logging.info")
    def test_logging_success(self, mock_info):
        """
        Checking that successful filtering is being logged.
        """
        input_fastq = "test_data/input_valid.fastq"
        output_fastq = "test_data/output_success.fastq"
        filter_fastq(input_fastq, output_fastq, (40, 60), (100, 500), 30)
        mock_info.assert_called_with(f"Все последовательности были отфильтрованы и сохранены в: '{output_fastq}'.")
        os.remove(output_fastq)
