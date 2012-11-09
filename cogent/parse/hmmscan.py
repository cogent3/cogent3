#!/usr/bin/env python

"""Parse output from hmmer3"""

def MinimalTbloutParser(lines):
    """Parse the .tblout output from hmmer3"""
    for line in lines:
        if line.startswith('#'):
            continue
        else:
            # lines are whitespace delimited, and there are 18 fields
            yield line.strip().split(None, 18)

headers_and_types = [
        ("target name", str),
        ("target accession", str), # renamed from accession due to ambiguity
        ("query name",str),
        ("query accession",str), # renamed from accession due to ambiguity
        ("full seq E-value",float), # renamed from E-value due to ambiguity
        ("full seq score",float), # renamed from score due to ambiguity
        ("full seq bias",float), # renamed from bias due to ambiguity
        ("best 1 domain E-value",float), # renamed from E-value due to ambiguity
        ("best 1 domain score",float), # renamed from score due to ambiguity
        ("best 1 domain bias",float), # renamed from bias due to ambiguity
        ("exp",float),
        ("reg",int),
        ("clu",int),
        ("ov",int),
        ("env",int),
        ("dom",int),
        ("rep",int),
        ("inc",int),
        ("description of target", str)]

def RichTbloutParser(lines):
    """Parse the .tblout output from hmmer3, return dicts per record"""
    for rec in MinimalTbloutParser(lines):
        yield dict([(h,t(r)) for (h,t),r in zip(headers_and_types, rec)])
