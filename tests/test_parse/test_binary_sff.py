#!/usr/bin/env python

import copy
import os
import tempfile
from unittest import TestCase, main

from cogent.parse.binary_sff import (
    seek_pad, parse_common_header, parse_read_header, parse_read_data,
    validate_common_header, parse_read, parse_binary_sff, UnsupportedSffError,
    write_pad, write_common_header, write_read_header, write_read_data,
    write_read, write_binary_sff, format_common_header, format_read_header,
    format_read_data, format_binary_sff, base36_encode, base36_decode,
    decode_location, decode_timestamp, decode_accession, decode_sff_filename,
    )

__author__ = "Kyle Bittinger"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Kyle Bittinger"]
__license__ = "GPL"
__version__ = "1.5.0"
__maintainer__ = "Kyle Bittinger"
__email__ = "kylebittinger@gmail.com"
__status__ = "Production"


TEST_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
SFF_FP = os.path.join(TEST_DIR, 'data', 'F6AVWTA01.sff')


class WritingFunctionTests(TestCase):
    def setUp(self):
        self.output_file = tempfile.TemporaryFile()

    def test_write_pad(self):
        self.output_file.write('\x01\x02\x03\x04')
        write_pad(self.output_file)
        self.output_file.seek(0)
        buff = self.output_file.read()
        self.assertEqual(buff, '\x01\x02\x03\x04\x00\x00\x00\x00')

    def test_write_common_header(self):
        write_common_header(self.output_file, COMMON_HEADER)

        file_pos = self.output_file.tell()
        self.assertTrue(file_pos % 8 == 0)

        self.output_file.seek(0)
        observed = parse_common_header(self.output_file)
        self.assertEqual(observed, COMMON_HEADER)

        file_pos = self.output_file.tell()
        self.assertTrue(file_pos % 8 == 0)

    def test_write_read_header(self):
        write_read_header(self.output_file, READ_HEADER)

        file_pos = self.output_file.tell()
        self.assertTrue(file_pos % 8 == 0)

        self.output_file.seek(0)
        observed = parse_read_header(self.output_file)
        self.assertEqual(observed, READ_HEADER)

        file_pos = self.output_file.tell()
        self.assertTrue(file_pos % 8 == 0)

    def test_write_read_data(self):
        write_read_data(self.output_file, READ_DATA)

        file_pos = self.output_file.tell()
        self.assertTrue(file_pos % 8 == 0)

        self.output_file.seek(0)
        num_flows = len(READ_DATA['flowgram_values'])
        num_bases = len(READ_DATA['Bases'])
        observed = parse_read_data(self.output_file, num_bases, num_flows)
        self.assertEqual(observed, READ_DATA)

        file_pos = self.output_file.tell()
        self.assertTrue(file_pos % 8 == 0)

    def test_write_read(self):
        read = READ_HEADER.copy()
        read.update(READ_DATA)
        write_read(self.output_file, read)

        file_pos = self.output_file.tell()
        self.assertTrue(file_pos % 8 == 0)

        self.output_file.seek(0)
        num_flows = len(read['flowgram_values'])
        observed = parse_read(self.output_file)
        self.assertEqual(observed, read)

        file_pos = self.output_file.tell()
        self.assertTrue(file_pos % 8 == 0)

    def test_write_binary_sff(self):
        read = READ_HEADER.copy()
        read.update(READ_DATA)

        header = COMMON_HEADER.copy()
        header['number_of_reads'] = 1

        write_binary_sff(self.output_file, header, [read])

        file_pos = self.output_file.tell()
        self.assertTrue(file_pos % 8 == 0)

        self.output_file.seek(0)
        observed_header, observed_reads = parse_binary_sff(
            self.output_file, native_flowgram_values=True)
        observed_reads = list(observed_reads)
        self.assertEqual(observed_header, header)
        self.assertEqual(observed_reads[0], read)
        self.assertEqual(len(observed_reads), 1)

        file_pos = self.output_file.tell()
        self.assertTrue(file_pos % 8 == 0)


class ParsingFunctionTests(TestCase):
    def setUp(self):
        self.sff_file = open(SFF_FP)

    def test_seek_pad(self):
        f = self.sff_file
        f.seek(8)
        seek_pad(f)
        self.assertEqual(f.tell(), 8)
        f.seek(9)
        seek_pad(f)
        self.assertEqual(f.tell(), 16)
        f.seek(10)
        seek_pad(f)
        self.assertEqual(f.tell(), 16)
        f.seek(15)
        seek_pad(f)
        self.assertEqual(f.tell(), 16)
        f.seek(16)
        seek_pad(f)
        self.assertEqual(f.tell(), 16)
        f.seek(17)
        seek_pad(f)
        self.assertEqual(f.tell(), 24)

    def test_parse_common_header(self):
        observed = parse_common_header(self.sff_file)
        self.assertEqual(observed, COMMON_HEADER)

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
        self.assertEqual(observed, READ_HEADER)

    def test_parse_read_data(self):
        self.sff_file.seek(440 + 32)
        observed = parse_read_data(self.sff_file, 271, 400)
        self.assertEqual(observed, READ_DATA)

    def test_parse_read(self):
        self.sff_file.seek(440)
        observed = parse_read(self.sff_file, 400)
        expected = dict(READ_HEADER.items() + READ_DATA.items())
        self.assertEqual(observed, expected)

    def test_parse_sff(self):
        header, reads = parse_binary_sff(self.sff_file)
        self.assertEqual(header, COMMON_HEADER)
        counter = 0
        for read in reads:
            self.assertEqual(
                len(read['flowgram_values']), header['number_of_flows_per_read'])
            counter += 1
        self.assertEqual(counter, 20)


class FormattingFunctionTests(TestCase):
    def setUp(self):
        self.output_file = tempfile.TemporaryFile()

    def test_format_common_header(self):
        self.assertEqual(
            format_common_header(COMMON_HEADER), COMMON_HEADER_TXT)

    def test_format_read_header(self):
        self.assertEqual(
            format_read_header(READ_HEADER), READ_HEADER_TXT)

    def test_format_read_header(self):
        self.assertEqual(
            format_read_data(READ_DATA, READ_HEADER), READ_DATA_TXT)

    def test_format_binary_sff(self):
        output_buffer = format_binary_sff(open(SFF_FP))
        output_buffer.seek(0)
        expected = COMMON_HEADER_TXT + READ_HEADER_TXT + READ_DATA_TXT
        observed = output_buffer.read(len(expected))
        self.assertEqual(observed, expected)


class Base36Tests(TestCase):
    def test_base36_encode(self):
        self.assertEqual(base36_encode(2), 'C')
        self.assertEqual(base36_encode(37), 'BB')

    def test_base36_decode(self):
        self.assertEqual(base36_decode('C'), 2)
        self.assertEqual(base36_decode('BB'), 37)

    def test_decode_location(self):
        self.assertEqual(decode_location('C'), (0, 2))

    def test_decode_timestamp(self):
        self.assertEqual(decode_timestamp('C3U5GW'), (2004, 9, 22, 16, 59, 10))
        self.assertEqual(decode_timestamp('GA202I'), (2010, 1, 22, 13, 28, 56))

    def test_decode_accession(self):
        self.assertEqual(
            decode_accession('GA202I001ER3QL'),
            ((2010, 1, 22, 13, 28, 56), '0', 1, (1843, 859)))

    def test_decode_sff_filename(self):
        self.assertEqual(
            decode_sff_filename('F6AVWTA01.sff'),
            ((2009, 11, 25, 14, 30, 19), 'A', 1))


COMMON_HEADER = {
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

COMMON_HEADER_TXT = """\
Common Header:
  Magic Number:  0x2E736666
  Version:       0001
  Index Offset:  33464
  Index Length:  900
  # of Reads:    20
  Header Length: 440
  Key Length:    4
  # of Flows:    400
  Flowgram Code: 1
  Flow Chars:    TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG
  Key Sequence:  TCAG
"""

READ_HEADER = {
    'name_length': 14,
    'Name': 'GA202I001ER3QL',
    'clip_adapter_left': 0,
    'read_header_length': 32,
    'clip_adapter_right': 0,
    'number_of_bases': 271,
    'clip_qual_left': 5,
    'clip_qual_right': 271,
    }

READ_HEADER_TXT = """
>GA202I001ER3QL
  Run Prefix:   R_2010_01_22_13_28_56_
  Region #:     1
  XY Location:  1843_0859

  Read Header Len:  32
  Name Length:      14
  # of Bases:       271
  Clip Qual Left:   5
  Clip Qual Right:  271
  Clip Adap Left:   0
  Clip Adap Right:  0
"""

READ_DATA = {
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

READ_DATA_TXT = """
Flowgram:	1.01	0.00	0.98	0.03	0.00	1.04	0.02	0.95	0.01	0.00	0.97	0.03	0.00	1.10	0.02	1.02	1.02	1.10	0.02	0.99	1.01	0.00	1.95	0.05	1.02	0.00	0.05	0.96	0.07	0.00	0.95	0.07	1.01	0.00	0.08	0.98	0.09	0.00	1.90	0.09	2.01	0.00	1.94	1.01	1.07	1.04	0.12	1.98	0.13	1.04	0.02	1.05	2.95	0.07	0.04	1.97	0.10	1.01	1.95	0.98	1.01	0.03	0.10	1.00	1.02	0.00	1.00	0.07	1.01	0.00	0.96	0.08	0.11	1.02	0.12	1.02	2.03	0.09	1.96	0.08	0.13	2.06	0.13	0.06	1.03	0.10	0.04	1.03	1.02	0.03	0.07	4.79	0.09	1.02	2.02	0.10	1.98	0.06	1.95	0.09	1.02	0.00	1.00	0.05	1.00	0.02	1.03	0.08	0.08	1.00	0.06	1.02	0.07	2.00	3.88	0.10	0.97	1.00	0.08	0.05	1.00	0.12	1.97	0.07	0.13	1.03	0.08	0.07	1.04	0.10	1.01	1.04	0.12	2.01	0.12	0.99	0.08	0.99	1.06	0.13	1.03	1.02	0.08	2.02	1.08	0.09	0.13	2.93	0.07	0.04	2.03	1.03	2.02	1.07	3.76	1.03	0.08	0.11	1.88	0.08	0.99	1.01	1.04	0.08	0.92	1.01	0.12	0.04	0.92	0.11	1.01	0.07	0.96	2.02	0.08	0.12	0.93	0.11	0.11	2.02	0.07	1.95	1.01	1.02	0.06	0.00	1.01	0.07	0.07	1.06	0.02	0.06	1.07	0.04	4.04	0.12	0.06	1.04	0.08	0.10	0.98	0.02	1.05	1.10	1.00	0.08	0.95	0.03	1.05	1.02	2.08	2.01	0.13	1.95	0.14	0.00	0.99	0.86	2.02	0.09	3.01	2.06	0.08	0.08	0.85	0.06	1.01	0.06	0.09	1.03	0.08	0.09	0.96	0.04	0.07	1.02	1.11	0.00	0.08	0.93	0.07	1.94	1.11	0.05	0.10	0.95	0.05	0.10	1.04	0.02	0.06	0.98	1.03	0.00	0.11	0.99	0.15	1.92	1.10	0.05	0.98	0.08	0.91	0.08	0.10	0.92	0.05	0.10	1.02	0.08	0.07	1.05	0.15	1.02	0.07	0.09	1.00	0.02	0.03	1.02	0.06	0.09	2.03	0.06	0.14	1.07	0.12	0.08	1.07	0.01	1.03	0.13	2.02	0.02	0.06	1.08	1.03	0.99	0.11	0.02	2.01	2.07	0.14	0.08	0.94	0.04	0.95	0.09	1.95	0.13	1.93	0.09	3.06	0.13	1.00	0.11	0.06	0.75	0.13	0.91	0.12	2.05	0.07	2.03	0.10	0.03	1.07	0.17	1.11	0.12	0.04	1.05	1.06	0.07	2.08	0.05	0.09	2.02	0.08	1.08	0.06	0.84	0.16	1.03	1.08	0.92	0.16	0.93	0.08	0.95	0.94	2.07	0.17	0.10	1.03	0.03	0.00	1.04	0.00	2.02	2.17	0.16	0.12	1.97	0.04	0.90	0.15	0.17	1.08	0.98	1.25	1.04	0.88	0.14	0.15	0.99	1.87	1.06	1.09	0.12	1.00	0.11	0.81	0.08	0.11	0.92	3.04	1.12	1.07	0.02	0.11	0.94	0.07	0.06	0.86	0.97	0.19	0.03	2.25	2.06
Flow Indexes:	1	3	6	8	11	14	16	17	18	20	21	23	23	25	28	31	33	36	39	39	41	41	43	43	44	45	46	48	48	50	52	53	53	53	56	56	58	59	59	60	61	64	65	67	69	71	74	76	77	77	79	79	82	82	85	88	89	92	92	92	92	92	94	95	95	97	97	99	99	101	103	105	107	110	112	114	114	115	115	115	115	117	118	121	123	123	126	129	131	132	134	134	136	138	139	141	142	144	144	145	148	148	148	151	151	152	153	153	154	155	155	155	155	156	159	159	161	162	163	165	166	169	171	173	174	174	177	180	180	182	182	183	184	187	190	193	195	195	195	195	198	201	203	204	205	207	209	210	211	211	212	212	214	214	217	218	219	219	221	221	221	222	222	225	227	230	233	236	237	240	242	242	243	246	249	252	253	256	258	258	259	261	263	266	269	272	274	277	280	283	283	286	289	291	293	293	296	297	298	301	301	302	302	305	307	309	309	311	311	313	313	313	315	318	320	322	322	324	324	327	329	332	333	335	335	338	338	340	342	344	345	346	348	350	351	352	352	355	358	360	360	361	361	364	364	366	369	370	371	372	373	376	377	377	378	379	381	383	386	387	387	387	388	389	392	395	396	399	399	400	400
Bases:	tcagCAGTAGTCCTGCTGCCTTCCGTAGGAGTTTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCGTCTACCGGAATACCTTTAATCATGTGAACATGTGAACTCATGATGCCATCTTGTATTAATCTTCCTTTCAGAAGGCTGTCCAAGAGTAGACGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGG
Quality Scores:	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	40	37	37	37	37	37	37	37	37	34	34	34	34	34	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	38	32	32	32	32	38	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	37	38	38	38	38	40	40	40	38	38	38	38	38	38	38	40	38	38	38	38	38	38	37	38	38	36	37	37	36	33	28	28	31	31	31	31	31	31	31	31	31	31	31	32	32	31	30	30	25	25	25	25
"""

if __name__ == '__main__':
    main()
