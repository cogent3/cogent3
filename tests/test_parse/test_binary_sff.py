import copy
import os
from unittest import TestCase, main

from cogent.parse.binary_sff import (
    pad, parse_common_header, parse_read_header, parse_read_data,
    validate_common_header, parse_read, parse_binary_sff, UnsupportedSffError,
    )

TEST_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
SFF_FP = os.path.join(TEST_DIR, 'data', 'F6AVWTA01.sff')


class FunctionTests(TestCase):
    def setUp(self):
        self.sff_file = open(SFF_FP)

    def test_pad(self):
        f = self.sff_file
        f.seek(8)
        pad(f)
        self.assertEqual(f.tell(), 8)
        f.seek(9)
        pad(f)
        self.assertEqual(f.tell(), 16)
        f.seek(10)
        pad(f)
        self.assertEqual(f.tell(), 16)
        f.seek(15)
        pad(f)
        self.assertEqual(f.tell(), 16)
        f.seek(16)
        pad(f)
        self.assertEqual(f.tell(), 16)
        f.seek(17)
        pad(f)
        self.assertEqual(f.tell(), 24)

    def test_parse_common_header(self):
        observed = parse_common_header(self.sff_file)
        self.assertEqual(observed, expected_common_header)

    def test_validate_common_header(self):
        header = {
            'magic_number': 779314790,
            'version': 1,
            'flowgram_format_code': 1,
            'index_offset': 0,
            'index_length': 0,
            'number_of_reads': 0,
            'header_length': 0,
            'key_length': 0,
            'number_of_flows_per_read': 0,
            'flow_chars': 'A',
            'key_sequence': 'A',
            }
        self.assertEqual(validate_common_header(header), None)
        header['version'] = 2
        self.assertRaises(UnsupportedSffError, validate_common_header, header)

    def test_parse_read_header(self):
        self.sff_file.seek(440)
        observed = parse_read_header(self.sff_file)
        self.assertEqual(observed, expected_read_header)

    def test_parse_read_data(self):
        self.sff_file.seek(440 + 32)
        observed = parse_read_data(self.sff_file, 271, 400)
        self.assertEqual(observed, expected_read_data)

    def test_parse_read(self):
        self.sff_file.seek(440)
        observed = parse_read(self.sff_file, 400)
        expected = dict(expected_read_header.items() + expected_read_data.items())
        self.assertEqual(observed, expected)

    def test_parse_sff(self):
        header, reads = parse_binary_sff(self.sff_file)
        self.assertEqual(header, expected_common_header)
        counter = 0
        for read in reads:
            self.assertEqual(
                len(read['flowgram_values']), header['number_of_flows_per_read'])
            counter += 1
        self.assertEqual(counter, 20)


expected_common_header = {
    'header_length': 440,
    'flowgram_format_code': 1,
    'index_length': 900,
    'magic_number': 779314790,
    'number_of_flows_per_read': 400,
    'version': 1,
    'flow_chars': 100 * 'TACG',
    'key_length': 4,
    'key_sequence': 'TCAG',
    'number_of_reads': 20,
    'index_offset': 33464,
    }

expected_read_header = {
    'name_length': 14,
    'Name': 'GA202I001ER3QL',
    'clip_adapter_left': 0,
    'read_header_length': 32,
    'clip_adapter_right': 0,
    'number_of_bases': 271,
    'clip_qual_left': 5,
    'clip_qual_right': 271,
    }

expected_read_data = {
    'flow_index_per_base': (
        1, 2, 3, 2, 3, 3, 2, 1, 1, 2, 1, 2, 0, 2, 3, 3, 2, 3, 3, 0, 2, 0, 2, 0,
        1, 1, 1, 2, 0, 2, 2, 1, 0, 0, 3, 0, 2, 1, 0, 1, 1, 3, 1, 2, 2, 2, 3, 2,
        1, 0, 2, 0, 3, 0, 3, 3, 1, 3, 0, 0, 0, 0, 2, 1, 0, 2, 0, 2, 0, 2, 2, 2,
        2, 3, 2, 2, 0, 1, 0, 0, 0, 2, 1, 3, 2, 0, 3, 3, 2, 1, 2, 0, 2, 2, 1, 2,
        1, 2, 0, 1, 3, 0, 0, 3, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 3, 0, 2, 1, 1, 2,
        1, 3, 2, 2, 1, 0, 3, 3, 0, 2, 0, 1, 1, 3, 3, 3, 2, 0, 0, 0, 3, 3, 2, 1,
        1, 2, 2, 1, 1, 0, 1, 0, 2, 0, 3, 1, 1, 0, 2, 0, 0, 1, 0, 3, 2, 3, 3, 3,
        1, 3, 2, 0, 1, 3, 3, 3, 1, 3, 2, 0, 1, 2, 2, 3, 3, 3, 2, 3, 3, 3, 0, 3,
        3, 2, 2, 0, 3, 1, 1, 3, 0, 1, 0, 3, 2, 2, 0, 2, 0, 2, 0, 0, 2, 3, 2, 2,
        0, 2, 0, 3, 2, 3, 1, 2, 0, 3, 0, 2, 2, 2, 1, 1, 2, 2, 1, 1, 0, 3, 3, 2,
        0, 1, 0, 3, 0, 2, 3, 1, 1, 1, 1, 3, 1, 0, 1, 1, 2, 2, 3, 1, 0, 0, 1, 1,
        3, 3, 1, 3, 0, 1, 0),
    'flowgram_values': (
        101, 0, 98, 3, 0, 104, 2, 95, 1, 0, 97, 3, 0, 110, 2, 102, 102, 110, 2,
        99, 101, 0, 195, 5, 102, 0, 5, 96, 7, 0, 95, 7, 101, 0, 8, 98, 9, 0,
        190, 9, 201, 0, 194, 101, 107, 104, 12, 198, 13, 104, 2, 105, 295, 7,
        4, 197, 10, 101, 195, 98, 101, 3, 10, 100, 102, 0, 100, 7, 101, 0, 96,
        8, 11, 102, 12, 102, 203, 9, 196, 8, 13, 206, 13, 6, 103, 10, 4, 103,
        102, 3, 7, 479, 9, 102, 202, 10, 198, 6, 195, 9, 102, 0, 100, 5, 100,
        2, 103, 8, 8, 100, 6, 102, 7, 200, 388, 10, 97, 100, 8, 5, 100, 12, 197,
        7, 13, 103, 8, 7, 104, 10, 101, 104, 12, 201, 12, 99, 8, 99, 106, 13,
        103, 102, 8, 202, 108, 9, 13, 293, 7, 4, 203, 103, 202, 107, 376, 103,
        8, 11, 188, 8, 99, 101, 104, 8, 92, 101, 12, 4, 92, 11, 101, 7, 96, 202,
        8, 12, 93, 11, 11, 202, 7, 195, 101, 102, 6, 0, 101, 7, 7, 106, 2, 6,
        107, 4, 404, 12, 6, 104, 8, 10, 98, 2, 105, 110, 100, 8, 95, 3, 105,
        102, 208, 201, 13, 195, 14, 0, 99, 86, 202, 9, 301, 206, 8, 8, 85, 6,
        101, 6, 9, 103, 8, 9, 96, 4, 7, 102, 111, 0, 8, 93, 7, 194, 111, 5, 10,
        95, 5, 10, 104, 2, 6, 98, 103, 0, 11, 99, 15, 192, 110, 5, 98, 8, 91, 8,
        10, 92, 5, 10, 102, 8, 7, 105, 15, 102, 7, 9, 100, 2, 3, 102, 6, 9, 203,
        6, 14, 107, 12, 8, 107, 1, 103, 13, 202, 2, 6, 108, 103, 99, 11, 2, 201,
        207, 14, 8, 94, 4, 95, 9, 195, 13, 193, 9, 306, 13, 100, 11, 6, 75, 13,
        91, 12, 205, 7, 203, 10, 3, 107, 17, 111, 12, 4, 105, 106, 7, 208, 5, 9,
        202, 8, 108, 6, 84, 16, 103, 108, 92, 16, 93, 8, 95, 94, 207, 17, 10,
        103, 3, 0, 104, 0, 202, 217, 16, 12, 197, 4, 90, 15, 17, 108, 98, 125,
        104, 88, 14, 15, 99, 187, 106, 109, 12, 100, 11, 81, 8, 11, 92, 304,
        112, 107, 2, 11, 94, 7, 6, 86, 97, 19, 3, 225, 206),
    'Bases': (
        'TCAGCAGTAGTCCTGCTGCCTTCCGTAGGAGTTTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCT'
        'CTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCC'
        'ATCGTCTACCGGAATACCTTTAATCATGTGAACATGTGAACTCATGATGCCATCTTGTATTAATCTTCCT'
        'TTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGG'),
    'quality_scores': (
        37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37,
        37, 37, 37, 37, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
        40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 37, 37, 37, 37, 37,
        37, 37, 37, 34, 34, 34, 34, 34, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37,
        37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37,
        37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37,
        38, 32, 32, 32, 32, 38, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37,
        37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37,
        37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37,
        37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37,
        37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37,
        37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37,
        37, 37, 37, 37, 38, 38, 38, 38, 40, 40, 40, 38, 38, 38, 38, 38, 38, 38,
        40, 38, 38, 38, 38, 38, 38, 37, 38, 38, 36, 37, 37, 36, 33, 28, 28, 31,
        31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 32, 32, 31, 30, 30, 25, 25, 25,
        25),
    }

if __name__ == '__main__':
    main()
