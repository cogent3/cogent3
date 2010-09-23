#!/usr/bin/env python
"""Parser for 454 Flowgram files in native binary format."""

__author__ = 'Kyle Bittinger'
__copyright__ = 'Copyright 2010, The Cogent Project'
__license__ = 'GPL'
__version__ = '1.5.0.dev'
__credits__ = ['Kyle Bittinger']
__maintainer__ = 'Kyle Bittinger'
__email__ = 'kylebittinger@gmail.com'
__status__ = 'Prototype'

import struct

# Inspired by, but not derived from, several other implementations:
# * BioPython (biopython.org)
# * sff_extract (www.melogen.upv.es/sff_extract)
# * Mothur (mothur.org)


class NamedStruct(struct.Struct):
    """Enhanced Struct class that associates names with each item in the struct.
    """

    def __init__(self, format, keys):
        """Create a new NamedStruct with a list of keys.
        """
        self.keys = keys
        super(NamedStruct, self).__init__(format)

    def read_from(self, file):
        """Read the struct from a file object and return the values as a dict.
        """
        buff = file.read(self.size)
        return self.unpack(buff)

    def pack(self, dict_of_vals):
        vals = [dict_of_vals[k] for k in self.keys]
        return super(NamedStruct, self).pack(*vals)

    def unpack(self, buffer):
        vals = super(NamedStruct, self).unpack(buffer)
        return dict(zip(self.keys, vals))


def seek_pad(file, unit=8):
    """Set a file's position to the next multiple of a given number.
    """
    position = file.tell()
    rem = position % unit
    if rem != 0:
        padding = unit - rem
        file.seek(padding, 1)        


def write_pad(file, unit=8):
    """Write zeros until the file's position is a multiple of the given number.
    """
    position = file.tell()
    rem = position % unit
    if rem != 0:
        num_bytes = unit - rem
        padding_bytes = '\x00' * num_bytes
        file.write(padding_bytes)


def parse_common_header(sff_file):
    """Parse a Common Header section from a binary SFF file.
    
    Keys in the resulting dict are identical to those defined in the
    Roche documentation.

    As a side effect, sets the position of the file object to the end
    of the Common Header section.
    """
    h = common_header_struct.read_from(sff_file)
    h['flow_chars'] = sff_file.read(h['number_of_flows_per_read'])
    h['key_sequence'] = sff_file.read(h['key_length'])
    seek_pad(sff_file)
    return h


def write_common_header(sff_file, header):
    """Write a common header section to a binary SFF file.
    """
    header_bytes = common_header_struct.pack(header)
    sff_file.write(header_bytes)
    sff_file.write(header['flow_chars'])
    sff_file.write(header['key_sequence'])
    write_pad(sff_file)


common_header_struct = NamedStruct('>IIQIIHHHB', [
    'magic_number',
    'version',
    'index_offset',
    'index_length',
    'number_of_reads',
    'header_length',
    'key_length',
    'number_of_flows_per_read',
    'flowgram_format_code',
    ])


class UnsupportedSffError(Exception):
    pass


def validate_common_header(header):
    """Validate the Common Header section of a binary SFF file.

    Raises an UnsupportedSffError if the header is not supported.
    """
    supported_values = {
        'magic_number': 0x2E736666,
        'version': 1,
        'flowgram_format_code': 1,
        }
    for attr_name, expected_value in supported_values.items():
        observed_value = header[attr_name]
        if observed_value != expected_value:
            raise UnsupportedSffError(
                '%s not supported. (Expected %s, observed %s)' % (
                    attr_name, expected_value, observed_value))


def parse_read_header(sff_file):
    """Parse a Read Header section from a binary SFF file.
    
    Keys in the resulting dict are identical to those defined in the
    Roche documentation.

    As a side effect, sets the position of the file object to the end
    of the Read Header section.
    """
    data = read_header_struct.read_from(sff_file)
    data['Name'] = sff_file.read(data['name_length'])
    seek_pad(sff_file)
    return data


def write_read_header(sff_file, header):
    """Write a read header section to a binary SFF file.
    """
    header_bytes = read_header_struct.pack(header)
    sff_file.write(header_bytes)
    sff_file.write(header['Name'])
    write_pad(sff_file)


read_header_struct = NamedStruct('>HHIHHHH', [
    'read_header_length',
    'name_length',
    'number_of_bases',
    'clip_qual_left', 
    'clip_qual_right',
    'clip_adapter_left',
    'clip_adapter_right',
    ])


def parse_read_data(sff_file, number_of_bases, number_of_flows=400):
    """Parse a Read Data section from a binary SFF file.
    
    Keys in the resulting dict are identical to those defined in the
    Roche documentation.

    As a side effect, sets the position of the file object to the end
    of the Read Header section.
    """
    data = {}
    flow_fmt = '>' + ('H' * number_of_flows)
    base_fmt = '>' + ('B' * number_of_bases)
    flow_fmt_size = struct.calcsize(flow_fmt)
    base_fmt_size = struct.calcsize(base_fmt)

    buff = sff_file.read(flow_fmt_size)
    data['flowgram_values'] = struct.unpack(flow_fmt, buff)

    buff = sff_file.read(base_fmt_size)
    data['flow_index_per_base'] = struct.unpack(base_fmt, buff)

    data['Bases'] = sff_file.read(number_of_bases)

    buff = sff_file.read(base_fmt_size)
    data['quality_scores'] = struct.unpack(base_fmt, buff)

    seek_pad(sff_file)
    return data


def write_read_data(sff_file, data):
    """Write a read data section to a binary SFF file.
    """
    number_of_flows = len(data['flowgram_values'])
    number_of_bases = len(data['quality_scores'])
    flow_fmt = '>' + ('H' * number_of_flows)
    base_fmt = '>' + ('B' * number_of_bases)

    flow_bytes = struct.pack(flow_fmt, *data['flowgram_values'])
    sff_file.write(flow_bytes)

    index_bytes = struct.pack(base_fmt, *data['flow_index_per_base'])
    sff_file.write(index_bytes)
    
    sff_file.write(data['Bases'])

    qual_bytes = struct.pack(base_fmt, *data['quality_scores'])
    sff_file.write(qual_bytes)

    write_pad(sff_file)


def parse_read(sff_file, number_of_flows=400):
    """Parse a single read from a binary SFF file.
    
    Keys in the resulting dict are identical to those defined in the
    Roche documentation for the Read Header and Read Data sections.

    As a side effect, sets the position of the file object to the end
    of the Read Data section.
    """
    header_data = parse_read_header(sff_file)
    read_data = parse_read_data(
        sff_file, header_data['number_of_bases'], number_of_flows)
    read_data.update(header_data)
    return read_data


def write_read(sff_file, read):
    """Write a single read to a binary SFF file.
    """
    write_read_header(sff_file, read)
    write_read_data(sff_file, read)


def parse_binary_sff(sff_file, native_flowgram_values=False):
    """Parse a binary SFF file, returning the header and a sequence of reads.

    In the binary file, flowgram values are stored as integers, 100
    times larger than the normalized floating point value.  Because
    the conversion is relatively expensive, we allow the computation
    to be skipped if the keyword argument native_flowgram_values is
    True.
    """
    header = parse_common_header(sff_file)
    number_of_flows = header['number_of_flows_per_read']
    validate_common_header(header)
    def get_reads():
        for i in range(header['number_of_reads']):

            # Skip the index section
            if sff_file.tell() == header['index_offset']:
                sff_file.seek(header['index_length'], 1)

            read = parse_read(sff_file, number_of_flows)

            if not native_flowgram_values:
                read['flowgram_values'] = [x * 0.01 for x in read['flowgram_values']]

            yield read
    return header, get_reads()

def write_binary_sff(sff_file, header, reads):
    """Write a binary SFF file, using provided header and read dicts.
    """
    sff_file.seek(0)
    sff_file.truncate()
    write_common_header(sff_file, header)
    for read in reads:
        write_read(sff_file, read)
