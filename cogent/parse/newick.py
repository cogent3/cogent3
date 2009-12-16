#!/usr/bin/env python
"""Newick format with all features as per the specs at:
    http://evolution.genetics.washington.edu/phylip/newick_doc.html
    http://evolution.genetics.washington.edu/phylip/newicktree.html
ie:
    Unquoted label underscore munging
    Quoted labels
    Inner node labels
    Lengths
    [ ... ] Comments (discarded)
    Unlabeled tips
also:
    Double quotes can be used.
    Spaces and quote marks are OK inside unquoted labels.
"""

from cogent.parse.record import FileFormatError
import re
EOT = None

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Peter Maxwell", "Andrew Butterfield", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"

class TreeParseError(FileFormatError):
    pass


class _Tokeniser(object):
    """An iterable stream of Newick tokens from 'text'
    
    By default this is very forgiving of non-standard unquoted labels.
    Two options can change how unquoted labels are interpreted:
      To prohibit internal spaces and quotes set strict_labels=True.
      To disable conversion of '_' to ' ' set underscore_unmunge=False.

    NOTE: underscore_unmunging is part of the Newick standard, although it
    is often inconvenient for other purposes.
    """
    
    def __init__(self, text, strict_labels=False, underscore_unmunge=True):
        self.text = text
        self.posn = None
        self.strict_unquoted_labels = strict_labels
        self.underscore_unmunge = underscore_unmunge
        self._generator = self._tokens()
        
    def __iter__(self):
        return self
        
    def next(self):
        return self._generator.next()
        
    def _errmsg(self):
        msg = 'Unexpected "%s"' % self.token
        if self.posn:
            (line, column) = self.posn
            if line == 'EOT':
                sample = line
            else:
                sample = self.text.split('\n')[line-1][
                        max(column-10, 0):column+len(self.token or '')-1]
            if line > 1:
                msg += ' at line %s:%s "%s"' % (line, column, sample)
            else:
                msg += ' at char %s "%s"' % (column, sample)
        return msg
        
    def _parse_error(self, detail=""):
        return TreeParseError(self._errmsg() + '. ' + detail)
                            
    def _tokens(self):
        closing_quote_token = None
        column = 1
        line = 1
        text = None
        closing_quote_token = None
        in_comment = False
        for token in re.split("""([\\t ]+|\\n|''|""|[]['"(),:;])""", self.text)+[EOT]:
            label_complete = False
            token_consumed = True
            
            self.token = token
            self.posn = (line, column)
            column += len(token or '')
            
            if token == "":
                pass
            elif in_comment:
                if token is EOT:
                    raise self._parse_error('Ended with unclosed comment')
                if token == ']':
                    in_comment = False
            elif closing_quote_token:
                if token is EOT:
                    raise self._parse_error('Text ended inside quoted label')
                if token == '\n':
                    raise self._parse_error('Line ended inside quoted label')
                if token == closing_quote_token:
                    label_complete = True
                    closing_quote_token = None
                else:
                    if token == closing_quote_token*2:
                        token = token[0]
                    text += token
            elif token is EOT or token in '\n[():,;':
                if text:
                    text = text.strip()
                    if self.underscore_unmunge:
                        text = text.replace('_', ' ')
                    label_complete = True
                if token == '\n':
                    line += 1
                    column = 1
                elif token == '[':
                    in_comment = True
                else:
                    token_consumed = False
            elif text is not None:
                text += token
            elif token in ["''", '""']:
                label_complete = True
                text = ""
            elif token in ["'", '"']:
                closing_quote_token = token
                text = ""
            elif token.strip():
                text = token
                label_complete = self.strict_unquoted_labels
            
            if label_complete:
                self.token = text
                yield text
                text = None
                
            if not token_consumed:
                self.token = token
                yield token  

def _parse_subtrees(tokens, constructor, sentinals):
    """A list of tree nodes from the 'tokens' generator.
    An ordinary newick file _should_ give a list of length 1.
    Calls constructor(children, name, attributes)"""
    nodes = []
    children = name = expected_attribute = None
    attributes = {}
    for token in tokens:
        if expected_attribute is not None:
            (attr_name, attr_cast) = expected_attribute
            try:
                attributes[attr_name] = attr_cast(token)
            except ValueError:
                raise tokens._parse_error("Can't convert %s '%s'" % 
                        (attr_name, token))
            expected_attribute = None
        elif token == '(':
            if children is not None:
                raise tokens._parse_error(
                    "Two subtrees in one node, missing comma?")
            elif name or attributes:
                raise tokens._parse_error(
                    "Subtree must be first element of the node.")
            children = _parse_subtrees(tokens, constructor, sentinals=[')'])
        elif token == ':':
            if 'length' in attributes:
                raise tokens._parse_error("Already have a length.")
            expected_attribute = ('length', float)
        elif token in [')', ';', ',', EOT]:
            if not (name or attributes or children or nodes) and token == ')':
                pass #  "();" as empty tree, not root with one unnamed tip.
            else:
                nodes.append(constructor(children, name, attributes))
            children = name = expected_attribute = None
            attributes = {}
            if token in sentinals:
                break
            elif token == ',' and ')' in sentinals:
                pass
            else:
                raise tokens._parse_error("Was expecting to end with %s" % 
                    ' or '.join([repr(s) for s in sentinals]))
        else:
            if name is not None:
                raise tokens._parse_error("Already have a name '%s' for this node." % name)
            elif attributes:
                raise tokens._parse_error("Name should come before length.")
            name = token
    return nodes
                
def parse_string(text, constructor, underscore_unmunge=True):
    """Parses a Newick-format string, using specified constructor for tree.

    Note: underscore_unmunge, if True, replaces underscores with spaces in
    the data that's read in. This is part of the Newick format, but it is
    often useful to suppress this behavior.
    """
    if text.strip()[0] not in  ["(", ";", ""]:
         # otherwise "filename" is a valid (if small) tree
        raise TreeParseError('Tree must start with "(", not "%s"' % text[:10])
    tokens = iter(_Tokeniser(text,underscore_unmunge=underscore_unmunge))
    trees = _parse_subtrees(tokens, constructor, [';', EOT])
    assert len(trees) ==  1, len(trees)
    return trees[0]
