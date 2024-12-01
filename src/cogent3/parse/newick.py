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

import re

from cogent3.parse.record import FileFormatError

EOT = None


class TreeParseError(FileFormatError):
    pass


class _Tokeniser:
    """Supplies an iterable stream of Newick tokens from 'text'

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

    def error(self, detail=""):
        if self.token:
            msg = f'Unexpected "{self.token}" at '
        else:
            msg = "At "
        (line, column) = self.posn
        sample = self.text.split("\n")[line][:column]
        if column > 30:
            sample = "..." + sample[-20:]
        if line > 0:
            msg += f'line {line + 1}:{column} "{sample}"'
        else:
            msg += f'char {column} "{sample}"'
        return TreeParseError(msg + ". " + detail)

    def tokens(self):
        column = 0
        line = 0
        text = None
        closing_quote_token = None
        in_comment = False
        comment = []
        splits = re.split(r"""([\t ]+|\n|''|""|[]['"(),:;\[\]])""", self.text)
        for token in splits + [EOT]:
            label_complete = False
            token_consumed = True
            self.token = token
            column += len(token or "")
            self.posn = (line, column)

            if token == "":  # we have to allow None to pass this
                continue

            if in_comment:
                if token is EOT:
                    raise self.error("Ended with unclosed comment")
                if token == "]":
                    in_comment = False
                    text = "".join(comment)
                    comment = []
                else:
                    comment.append(token)
                    continue

            if closing_quote_token:
                if token is EOT:
                    raise self.error("Text ended inside quoted label")
                if token == "\n":
                    raise self.error("Line ended inside quoted label")
                if token == closing_quote_token:
                    label_complete = True
                    closing_quote_token = None
                else:
                    if token == closing_quote_token * 2:
                        token = token[0]
                    text += token
            elif token is EOT or token in "\n[]():,;":
                if text:
                    text = text.strip()
                    if self.underscore_unmunge and "_" in text:
                        text = text.replace("_", " ")
                    label_complete = True
                if token == "\n":
                    line += 1
                    column = 1
                elif token == "[":
                    in_comment = True
                    token_consumed = False
                elif token != "]":
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
                self.token = None
                yield text
                text = None

            if not token_consumed:
                self.token = token
                yield token


def parse_string(text, constructor, **kw):
    """Parses a Newick-format string, using specified constructor for tree.

    Calls constructor(children, name, attributes)

    Note: underscore_unmunge, if True, replaces underscores with spaces in
    the data that's read in. This is part of the Newick format, but it is
    often useful to suppress this behavior.
    """
    if "(" not in text and ";" not in text and text.strip():
        # otherwise "filename" is a valid (if small) tree
        raise TreeParseError(f'Not a Newick tree: "{text[:10]}"')
    sentinals = [";", EOT]
    stack = []
    nodes = []
    children = name = expected_attribute = None
    attributes = {}
    tokeniser = _Tokeniser(text, **kw)
    for token in tokeniser.tokens():
        if expected_attribute is not None:
            (attr_name, attr_cast) = expected_attribute
            try:
                attributes[attr_name] = attr_cast(token)
            except ValueError:
                raise tokeniser.error(f"Can't convert {attr_name!r} '{token}'")
            expected_attribute = None
        elif token == "(":
            if children is not None:
                raise tokeniser.error("Two subtrees in one node, missing comma?")
            if name or attributes:
                raise tokeniser.error("Subtree must be first element of the node.")
            stack.append((nodes, sentinals, attributes))
            (nodes, sentinals, attributes) = ([], [")"], {})
        elif token == ":":
            if "length" in attributes:
                raise tokeniser.error("Already have a length.")
            expected_attribute = ("length", float)
        elif token == "[":
            if "other" in attributes:
                raise tokeniser.error("Already have a 'other'.")
            expected_attribute = ("other", lambda x: x.split(","))
        elif token in [")", ";", ",", EOT]:
            nodes.append(constructor(children, name, attributes))
            children = name = expected_attribute = None
            attributes = {}
            if token in sentinals:
                if not stack:
                    break
                children = nodes
                (nodes, sentinals, attributes) = stack.pop()
            elif token != "," or ")" not in sentinals:
                raise tokeniser.error(
                    f'Was expecting to end with {" or ".join([repr(s) for s in sentinals])}',
                )
        elif name is None:
            name = token
        else:
            raise tokeniser.error(f"Already have a name {name!r} for this node.")
    assert not stack, stack
    assert len(nodes) == 1, len(nodes)
    return nodes[0]
