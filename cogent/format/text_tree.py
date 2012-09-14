#!/bin/env python
# file text_tree.py
"""Simple base text representation of phylo tree."""

__author__ = "Micah Hamady"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Micah Hamady"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Micah Hamady"
__email__ = "hamady@colorado.edu"
__status__ = "Prototype"

class TextNode(object):
    """ Helper to display text phyolo tree """
    def __init__(self, branch_length, name):
        """Set up defaults"""
        self.BranchLength = 0.0
        if branch_length:
            self.BranchLength = branch_length

        self.Name = name
        self.Printed = False
        self.NumChildren = 0
        self.NumChildrenPrinted = 0
        self.LastChild = None

    def display(self, max_dist, scale):
        """Display current node - should refactor this"""
        delimiter = "-"
        pdelimiter = "+"

        if self.Printed:
            delimiter = " "
            pdelimiter = "|"

        # compute number of children last child has
        last_ct = 0
        if self.LastChild:
            last_ct = self.LastChild.NumChildren

        # update values
        self.Printed = True
        self.NumChildrenPrinted += 1
        print_ct = self.NumChildren - last_ct

        if self.NumChildrenPrinted == (print_ct + 1): # or \
            delimiter = " "
            pdelimiter = "+"
        elif self.NumChildrenPrinted > (print_ct + 1):
            pdelimiter = " "
            delimiter = " "
        if (self.NumChildren == self.NumChildrenPrinted and self.NumChildren == print_ct): 
            delimiter = " "
            pdelimiter = "+"
 
        # check if leaf
        dout = ""
        if self.Name:
            dout = "@@%s" % self.Name 
            pdelimiter = ">"
            delimiter = "-"

        return (int(self.BranchLength) - 1) * delimiter + pdelimiter + dout

def process_nodes(all_nodes, cur_nodes, cur_node, parent):
    """ Recursively process nodes """
    # make current node
    pn = TextNode(cur_node.Length, cur_node.Name)

    # set current node as last child of parent (last one wins)
    if parent:
        parent.LastChild = pn

    # handle terminal node
    if not cur_node.Children:
        all_nodes.append(cur_nodes + [pn])

    # internal node
    else:
        cur_nodes.append(pn)
        pn.NumChildren = 0 
        for child in cur_node.Children:
            if child.Children:
                pn.NumChildren +=  len(child.tips())
            else:
                pn.NumChildren += 1 
        for child in cur_node.Children: 
            process_nodes(all_nodes, cur_nodes, child, pn)
        cur_nodes.pop()

def generate_nodes(tree, max_dist, scale):
    """Iterate over list of TextNodes to display """ 
    all_nodes = []
    cur_nodes = []
    
    # generate list of lists of TextNodes
    process_nodes(all_nodes, cur_nodes, tree, None)

    # process each list of TextNodes
    for node_list in all_nodes:

        # generate text string and node key 
        ticks, node_key =  ''.join([x.display(max_dist, scale) for x in node_list]).split("@@")
        
        # compute distances
        dist = sum([x.BranchLength for x in node_list])
        #scaled_dist = sum([x.ScaledBranchLength for x in node_list])
        scaled_dist = sum([x.BranchLength for x in node_list])
        branch_len = node_list[-1].BranchLength

        # yield each node 
        yield ticks, node_key, dist, scaled_dist, branch_len

