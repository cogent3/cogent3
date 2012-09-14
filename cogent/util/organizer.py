#!/usr/bin/env python
"""Provides Filter, Organizer, GroupList.

Filter is a dictionary of {field_name:[filter_functions]}.
Organizer is initiated with a list of Filters and it is called on a list of
objects. It organizes all objects according to the given criteria (Filters).
GroupList stores a group of objects and a separte hierarchy (i.e. list
of elements or functions it matched).
"""

__author__ = "Sandra Smit"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Sandra Smit"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Sandra Smit"
__email__ = "sandra.smit@colorado.edu"
__status__ = "Development"

class Filter(dict):
    """Dictionary of functions, i.e. selection criteria"""

    def __init__(self,Name,functions_dict=None):
        """Returns a new Filter object with given name and funtions"""
        if not functions_dict:
            functions_dict = {}
        super(Filter, self).__init__(functions_dict)
        self.Name = Name

    def __call__(self, item):
        """Returns True if the item satisfies all criteria"""
        for field in self:
            for function in self[field]:
                try:
                    if field is None:
                        if not function(item):
                            return False
                    else:
                        if not function(getattr(item, field)):
                            return False
                except AttributeError:
                    return False
        return True

class GroupList(list):
    """GroupList is a list of objects accompanied by a hierarchy.

    The data can be of any type.
    Groups is a list that stores the hierarchy of things the data matched.
    E.g. filters, classifiers.

    Example GroupList: data = [1,2,3,4], groups = ['numbers','<5']
            GroupList: data = [5,6,7,8], groups = ['numbers','>=5']
    """

    def __init__(self, data, Groups=None):
        """Returns a new GroupList with data and Groups"""
        super(GroupList, self).__init__(data)
        self.Groups = Groups or []


class Organizer(list):
    """Organizer puts given objects into dict according to given criteria"""

    def _find_first_match(self, item):
        """Returns the name of the first criterion that applies to the item"""
        for criterion in self:
            if criterion(item):
                return criterion.Name

    def __call__(self, data, existing_results=None):
        """Puts the sequences into the dictionary according to the criteria"""
        temp_results = {None:[]}   #dict to store intermediate results
        results = []        #the final results
        
        try:
            groups = data.Groups
        except AttributeError:
            groups = []

        #figure out what criteria we have, and make keys for them
        for criterion in self:
            temp_results[criterion.Name] = []

        #partition each datum according to which filter it matches first
        for item in data:
            temp_results[self._find_first_match(item)].append(item)
        
        return [GroupList(value, groups + [key]) \
            for key, value in temp_results.items() if value]

def regroup(data):
    """Regroups data into matching categories.
    
    data is a list of GroupLists (a list with hierarchy info associated)
    e.g. regroup([GroupList([1,2],['a']),GroupList([8,9],['b']),
            GroupList([3,4],['a'])]
    will give: [[1,2,3,4],[8,9]] 
    """
    result = {}
    for group_list in data:
        analyses = tuple(group_list.Groups)
        
        try:
            result[analyses].extend(group_list)
        except:
            result[analyses] = []
            result[analyses].extend(group_list)
    return  [GroupList(result[k],list(k)) for k,v in result.items()]
