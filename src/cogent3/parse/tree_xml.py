#!/usr/bin/env python
"""Parses a simple but verbose XML representation of a phylogenetic tree, with
elements <clade>, <name>, <param> and <value>.  XML attributes are not used,
so the syntax of parameter names is not restricted at all.

Newick
------
((a,b:3),c);

XML
---
<clade>
  <clade>
    <clade>
       <name>a</name>
    </clade>
    <clade>
       <name>b</name>
       <param><name>length</name><value>3.0</value></param>
    </clade>
  </clade>
  <clade>
     <name>c</name>
  </clade>
</clade>

Parameters are inherited by contained clades unless overridden.

"""

import xml.sax


__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"


class TreeHandler(xml.sax.ContentHandler):
    def __init__(self, tree_builder):
        self.build_edge = tree_builder

    def startDocument(self):
        self.stack = [({}, None, None)]
        self.data = {"clades": [], "params": {}}
        self.in_clade = False
        self.current = None

    def startElement(self, name, attrs):
        self.parent = self.data
        self.stack.append((self.data, self.in_clade, self.current))
        self.current = ""
        if name == "clade":
            self.data = {
                "params": self.data["params"].copy(),
                "clades": [],
                "name": None,
            }
            self.in_clade = True
        else:
            self.data = {}
            self.in_clade = False

    def characters(self, text):
        self.current += str(text)

    def endElement(self, name):
        getattr(self, f"process_{name}")(self.current, **self.data)
        (self.data, self.in_clade, self.current) = self.stack.pop()
        self.parent = self.stack[-1][0]

    def endDocument(self):
        pass

    def process_clade(self, text, name, params, clades):
        edge = self.build_edge(clades, name, params)
        self.parent["clades"].append(edge)

    def process_param(self, text, name, value):
        self.parent["params"][name] = value

    def process_name(self, text):
        self.parent["name"] = text.strip()

    def process_value(self, text):
        if text == "None":
            self.parent["value"] = None
        else:
            self.parent["value"] = float(text)


def parse_string(text, tree_builder):
    handler = TreeHandler(tree_builder)
    xml.sax.parseString(text.encode("utf8"), handler)
    trees = handler.data["clades"]
    assert len(trees) == 1, trees
    return trees[0]
