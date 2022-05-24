# /usr/bin/env python
"""Parsers for tree formats.

Implementation Notes

The algorithm used here is fairly general: should possibly make the code
generalizable to tree strings that use alternative delimiters and symbols.
However, I can't think of any cases where alternatives are used, so this is
left to future work.

Should possibly build a dict of {label:TreeNode} while parsing to make it
convenient to fill in additional data later, e.g. to fill in sequences from
their numeric labels in Newick format. Alternatively, maybe TreeNode should
get a buildIndex() method that performs the equivalent task.

As of 12/27/03, should be capable of parsing the ClustalW .dnd files without
difficulty.

"""
from cogent3.core.tree import PhyloNode
from cogent3.parse.record import RecordError


__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Rob Knight", "Catherine Lozupone", "Daniel McDonald"]
__license__ = "BSD-3"
__version__ = "2022.5.25a1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Development"

strip = str.strip
maketrans = str.maketrans

_dnd_token_str = "(:),;"
_dnd_tokens = dict.fromkeys(_dnd_token_str)
_dnd_tokens_and_spaces = _dnd_token_str + " \t\v\n"

remove_dnd_tokens = maketrans(_dnd_tokens_and_spaces, "-" * len(_dnd_tokens_and_spaces))


def safe_for_tree(s):
    """Makes string s safe for DndParser by removing significant chars."""
    return s.translate(remove_dnd_tokens)


def bad_dnd_tokens(s, is_valid_name):
    """Returns list of bad dnd tokens from s, using is_valid_name for names.

    Useful for finding trees with misformatted names that break parsing.
    """
    for t in DndTokenizer(s):
        if t in _dnd_tokens:
            continue
        # also OK if it's a number
        try:
            float(t)
            continue
        except:  # wasn't a number -- further tests
            pass
        if is_valid_name(t):
            continue
        # if we got here, nothing worked, so yield the current token
        yield t


def DndTokenizer(data):
    """Tokenizes data into a stream of punctuation, labels and lengths.

    Note: data should all be a single sequence, e.g. a single string.
    """
    in_quotes = False
    saved = []
    sa = saved.append
    for d in data:
        if d == "'":
            in_quotes = not (in_quotes)
        if d in _dnd_tokens and not in_quotes:
            curr = "".join(saved).strip()
            if curr:
                yield curr
            yield d
            saved = []
            sa = saved.append
        else:
            sa(d)


def DndParser(lines, constructor=PhyloNode, unescape_name=False):
    """Returns tree from the Clustal .dnd file format, and anything equivalent.

    Tree is made up of cogent3.base.tree.PhyloNode objects, with branch lengths
    (by default, although you can pass in an alternative constructor
    explicitly).
    """
    if isinstance(lines, str):
        data = lines
    else:
        data = "".join(lines)
    # skip arb comment stuff if present: start at first paren
    paren_index = data.find("(")
    data = data[paren_index:]
    left_count = data.count("(")
    right_count = data.count(")")
    if left_count != right_count:
        raise RecordError(
            f"Found {left_count} left parens but {right_count} right parens."
        )

    tokens = DndTokenizer(data)
    curr_node = None
    state = "PreColon"
    state1 = "PreClosed"
    last_token = None
    for t in tokens:
        if t == ":":  # expecting branch length
            state = "PostColon"
            # prevent state reset
            last_token = t
            continue
        if t == ")" and (last_token == "," or last_token == "("):  # node without name
            new_node = _new_child(curr_node, constructor)
            new_node.name = None
            curr_node = new_node.parent
            state1 = "PostClosed"
            last_token = t
            continue
        if t == ")":  # closing the current node
            curr_node = curr_node.parent
            state1 = "PostClosed"
            last_token = t
            continue
        if t == "(":  # opening a new node
            curr_node = _new_child(curr_node, constructor)
        elif t == ";":  # end of data
            last_token = t
            break
        # node without name
        elif t == "," and (last_token == "," or last_token == "("):
            new_node = _new_child(curr_node, constructor)
            new_node.name = None
            curr_node = new_node.parent
        elif t == ",":  # separator: next node adds to this node's parent
            curr_node = curr_node.parent
        elif state == "PreColon" and state1 == "PreClosed":  # data for the current node
            new_node = _new_child(curr_node, constructor)
            if unescape_name:
                if t.startswith("'") and t.endswith("'"):
                    while t.startswith("'") and t.endswith("'"):
                        t = t[1:-1]
                else:
                    if "_" in t:
                        t = t.replace("_", " ")
            new_node.name = t
            curr_node = new_node
        elif state == "PreColon" and state1 == "PostClosed":
            if unescape_name:
                while t.startswith("'") and t.endswith("'"):
                    t = t[1:-1]
            curr_node.name = t
        elif state == "PostColon":  # length data for the current node
            curr_node.length = float(t)
        else:  # can't think of a reason to get here
            raise RecordError(f"Incorrect PhyloNode state? {t}")
        state = "PreColon"  # get here for any non-colon token
        state1 = "PreClosed"
        last_token = t

    if curr_node is not None and curr_node.parent is not None:
        raise RecordError("Didn't get back to root of tree.")

    if curr_node is None:  # no data -- return empty node
        return constructor()
    return curr_node  # this should be the root of the tree


def _new_child(old_node, constructor):
    """Returns new_node which has old_node as its parent."""
    new_node = constructor()
    new_node.parent = old_node
    if old_node is not None:
        if id(new_node) not in list(map(id, old_node.children)):
            old_node.children.append(new_node)
    return new_node
