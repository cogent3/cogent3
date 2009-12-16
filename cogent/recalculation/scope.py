#!/usr/bin/env python
from __future__ import division
import logging
from cogent.util.dict_array import DictArrayTemplate
from cogent.maths.stats.distribution import chdtri

LOG = logging.getLogger('cogent')

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"

class ScopeError(KeyError):
    pass

class InvalidScopeError(ScopeError):
    """for scopes including an unknown value for a known dimension"""
    pass

class InvalidDimensionError(ScopeError):
    """for scopes including an unknown dimension"""
    pass

class IncompleteScopeError(ScopeError):
    """For underspecified scope when retrieving values"""
    pass


# Can be passed to _LeafDefn.interpretScopes()
class _ExistentialQualifier(object):
    def __init__(self, cats=None):
        self.cats = cats
    def __repr__(self):
        if self.cats is None:
            return self.__class__.__name__
        else:
            return '%s(%s)' % (self.__class__.__name__, self.cats)
    

class EACH(_ExistentialQualifier):
    independent = True

class ALL(_ExistentialQualifier):
    independent = False


def theOneItemIn(items):
    assert len(items) == 1, items
    return iter(items).next()

def _indexed(values):
    # This is the core of the redundancy elimination, used to group
    # identical calculations.
    # >>> _indexed({'a':1.0, 'b':2.0, 'c':3.0, 'd':1.0, 'e':1.0})
    # ([1.0, 2.0, 3.0], {'a':0, 'b':1, 'c':2, 'd':0, 'e':0})
    uniq = []
    index = {}
    values = values.items()
    values.sort()
    for (key, value) in values:
        if value in uniq:
            u = uniq.index(value)
        else:
            u = len(uniq)
            uniq.append(value)
        index[key] = u
    return uniq, index


def _fmtrow(width, values, maxwidth):
    if (len(dict([(id(v),1) for v in values])) == 1 and
            len(str(values[0])) > width):
        s = str(values[0]).replace('\n', ' ')
        if len(s) > maxwidth:
            s = s[:maxwidth-4] + '...'
    else:
        template = '%%%ss' % width
        s = ''.join([(template % (v,)).replace('\n', ' ')[:width]
                for v in values])
    return s


class Undefined(object):
    # Placeholder for a value that can't be calculated
    # because input 'name' has not been provided.
    def __init__(self, name):
        self.name = name
    def __repr__(self):
        return 'Undef(%s)' % self.name
    

def nullor(name, f, recycled=False):
    # If None, record as undefined.
    # If undefined, propagate error up.
    # Otherwise, do the calculation.
    def g(*args):
        undef = [x for x in args if isinstance(x, Undefined)]
        if undef:
            return undef[0]
        elif None in args:
            return Undefined(name)
        else:
            if recycled:
                args = (None,) + args
            return f(*args)
    return g

#  Level1:  D E F I N I T I O N S

#  Each ParamDefn supplied a .calc(args) method. Used to define the
#  calculation as a DAG of ParamDefns.

# A _Defn has two phases in its life: pre activation it just has .args,
# post activation (once it becomes part of a parameter controller) it
# holds a dynamic list of scope assignments.
# This means defn.makeParamController() can only be called once.

class _Defn(object):
    name = '?'
    default = None
    user_param = False
    
    def __init__(self):
        self.clients = []
        self.selection = {}
        self.assignments = {}
        self.activated = False
    
    def makeName(self, name, extra_label=None):
        if name is None:
            name = self.name
        if extra_label is not None:
            name += extra_label
        return name
    
    def getDefaultSetting(self):
        return None
    
    def addClient(self, client):
        assert not self.activated, self.name
        assert not self.assignments, self.assignments
        self.clients.append(client)
    
    def acrossDimension(self, dimension, cats):
        return [self.selectFromDimension(dimension, cat) for cat in cats]
    
    def selectFromDimension(self, dimension, cat):
        return SelectFromDimension(self, **{dimension:cat})
    
    def getRequiredScopes(self, arg_dimensions):
        # A list of scope dictionaries: [{dimension:value},] that this
        # Defn needs from an input Defn with `arg_dimensions`
        if not self.activated:
            assert not self.clients, self.clients
            raise RuntimeError('Value at "%s" step never used' % self.name)
        if self.assignments:
            result = []
            for scope_t in self.assignments:
                sel = {}
                sel.update(self.selection)
                for (d, c) in zip(self.valid_dimensions, scope_t):
                    if d in arg_dimensions:
                        sel[d] = c
                result.append(sel)
        else:
            result = [self.selection]
        return result
    
    def addScopes(self, scopes):
        assert not self.activated
        for scope in scopes:
            scope_t = [scope.get(d, 'all') for d in self.valid_dimensions]
            scope_t = tuple(scope_t)
            if scope_t not in self.assignments:
                self.assignments[scope_t] = self.getDefaultSetting()
    
    def outputOrdinalFor(self, scope):
        scope_t = tuple([scope[d] for d in self.valid_dimensions])
        return self.index[scope_t]
    
    def usedDimensions(self):
        used = []
        for (d, dim) in enumerate(self.valid_dimensions):
            seen = {}
            for (scope_t, i) in self.index.items():
                rest_of_scope = scope_t[:d]+scope_t[d+1:]
                if rest_of_scope in seen:
                    if i != seen[rest_of_scope]:
                        used.append(dim)
                        break
                else:
                    seen[rest_of_scope] = i
        return tuple(used) + self.internal_dimensions
    
    def _getPosnForScope(self, *args, **scope):
        scope = self.interpretPositionalScopeArgs(*args, **scope)
        posns = set()
        for scope_t in self.interpretScope(**scope):
            posns.add(self.index[scope_t])
        if len(posns) == 0:
            raise InvalidScopeError("no value for %s at %s" % (self.name, scope))
        if len(posns) > 1:
            raise IncompleteScopeError("%s distinct values of %s within %s" %
                    (len(posns), self.name, scope))
        return theOneItemIn(posns)
    
    def wrapValue(self, value):
        if isinstance(value, Undefined):
            raise ValueError('Input "%s" is not defined' % value.name)
        if getattr(self, 'array_template', None) is not None:
            value = self.array_template.wrap(value)
        return value
    
    def unwrapValue(self, value):
        if getattr(self, 'array_template', None) is not None:
            value = self.array_template.unwrap(value)
        return value
    
    def getCurrentValueForScope(self, *args, **scope):
        posn = self._getPosnForScope(*args, **scope)
        return self.wrapValue(self.values[posn])
    
    def getCurrentSettingForScope(self, *args, **scope):
        posn = self._getPosnForScope(*args, **scope)
        return self.uniq[posn]
    
    def interpretPositionalScopeArgs(self, *args, **scope):
        # Carefully turn scope args into scope kwargs
        assert len(args) <= len(self.valid_dimensions), args
        for (dimension, arg) in zip(self.valid_dimensions, args):
            assert dimension not in scope, dimension
            scope[dimension] = arg
        return scope
    
    def interpretScopes(self, independent=None, **kw):
        """A list of the scopes defined by the selecting keyword arguments.
        
        Keyword arguments should be of the form dimension=settings,
        where settings are a list of categories from that
        dimension, or an instance of EACH or ALL wrapping such a list.
        
        A missing list, None, or an uninstantiated ALL / EACH class
        is taken to mean the entire dimension.
        
        If 'independent' (which defaults to self.independent_by_default)
        is true then category lists not wrapped as an EACH or an ALL will
        be treated as an EACH, otherwise as an ALL.
        
        There will only be one scope in the resulting list unless at least
        one dimension is set to EACH."""
        
        if independent is None:
            independent = self.independent_by_default
        
        # interpretScopes is used for assigning, so should specify
        # the scope exactly
        for d in kw:
            if d not in self.valid_dimensions:
                raise InvalidDimensionError(d)
        
        # Initially ignore EACH, just get a full ungrouped set
        kw2 = {}
        independent_dimensions = []
        for (i,dimension) in enumerate(self.valid_dimensions):
            selection = kw.get(dimension, None)
            if selection in [EACH, ALL]:
                dimension_independent = selection.independent
                selection = None
            elif isinstance(selection, (EACH, ALL)):
                dimension_independent = selection.independent
                selection = selection.cats
            else:
                dimension_independent = independent
            if dimension_independent:
                independent_dimensions.append(i)
            if selection is not None:
                kw2[dimension] = selection
        all = self.interpretScope(**kw2)
        
        # Group independent scopes
        result = {}
        for scope_t in all:
            key = tuple([scope_t[i] for i in independent_dimensions])
            if key not in result:
                result[key] = set()
            result[key].add(scope_t)
        return result.values()
    
    def interpretScope(self, **kw):
        """A set of the scope-tuples that match the input dict like
        {dimension:[categories]}"""
        selector = []
        unused = {}
        valid_dimensions = list(self.valid_dimensions)
        for d in kw:
            if d not in valid_dimensions:
                continue
            if kw[d] is None: #i.e.: ALL
                continue
            if isinstance(kw[d], str):
                kw[d] = [kw[d]]
            assert type(kw[d]) in [tuple, list], (d, kw[d])
            assert len(kw[d]), kw[d]
            selector.append((valid_dimensions.index(d), d, kw[d]))
            unused[d] = kw[d][:]
        
        result = set()
        for scope_t in self.assignments:
            for (i, d, cs) in selector:
                if scope_t[i] not in cs:
                    break
            else:
                result.add(scope_t)
                for (i, d, cs) in selector:
                    if d in unused:
                        if scope_t[i] in unused[d]:
                            unused[d].remove(scope_t[i])
                        if not unused[d]:
                            del unused[d]
        
        if unused:
            # print unused, self.assignments.keys()
            raise InvalidScopeError(unused)
        
        return result
    
    def fillParValueDict(self, result, dimensions, cell_value_lookup):
        
        assert self.name not in result, self.name
        
        posns = [
                list(self.valid_dimensions).index(d)
                for d in dimensions
                if d in self.valid_dimensions]
        for (scope_t, i) in self.index.items():
            value = cell_value_lookup(self, i)
            value = self.wrapValue(value)
            scope = tuple([scope_t[i] for i in posns])
            
            (d,key) = (result, self.name)
            for key2 in scope:
                if key not in d: d[key] = {}
                (d, key) = (d[key], key2)
            
            if key in d and  value != d[key]:
                msg = 'Multiple values for %s' % self.name
                if scope:
                    msg += ' within scope %s' % '/'.join(scope)
                raise IncompleteScopeError(msg)
            d[key] = value
    
    def _update_from_assignments(self):
        (self.uniq, self.index) = _indexed(self.assignments)
    
    def _local_repr(self, col_width, max_width):
        body = []
        for (i, arg) in enumerate(self.args):
            row = []
            if isinstance(arg, SelectFromDimension):
                argname = arg.arg.name
                for nums in self.uniq:
                    num = arg.uniq[nums[i]]
                    row.append(theOneItemIn(num))
            else:
                argname = arg.name
                for nums in self.uniq:
                    row.append(nums[i])
            body.append((['', self.name][i==0], argname, row))
        
        return '\n'.join(
            ['%-10s%-10s%s' % (label1[:9], label2[:9],
                    _fmtrow(col_width+1, settings, max_width))
            for (label1, label2, settings) in body])
    
    def __repr__(self):
        return '%s(%s x %s)' % (self.__class__.__name__, self.name,
                len(getattr(self, 'cells', [])))
    

class SelectFromDimension(_Defn):
    """A special kind of Defn used to bridge from Defns where a particular
    dimension is just part of the scope rules to later Defns where each
    value has its own Defn, eg: edges of a tree"""
    
    name = 'select'
    #params = {}
    
    def __init__(self, arg, **kw):
        assert not arg.activated, arg.name
        _Defn.__init__(self)
        self.args = (arg,)
        self.arg = arg
        self.valid_dimensions = tuple([
            d for d in arg.valid_dimensions if d not in kw])
        self.selection = kw
        arg.addClient(self)
    
    def update(self):
        for scope_t in self.assignments:
            scope = dict(zip(self.valid_dimensions, scope_t))
            scope.update(self.selection)
            input_num = self.arg.outputOrdinalFor(scope)
            self.assignments[scope_t] = (input_num,)
        self._update_from_assignments()
        self.values = [self.arg.values[i] for (i,) in self.uniq]
    
    def makeCells(self, input_soup, variable=None):
        arg = input_soup[id(self.arg)]
        outputs = [arg[input_num] for (input_num,) in self.uniq]
        return ([], outputs)     
    

class _NonLeafDefn(_Defn):
    def __init__(self, *args, **kw):
        _Defn.__init__(self)
        valid_dimensions = []
        for arg in args:
            assert isinstance(arg, _Defn), type(arg)
            assert not arg.activated, arg.name
            for dimension in arg.valid_dimensions:
                if dimension not in valid_dimensions:
                    valid_dimensions.append(dimension)
            #print >>sys.stderr, arg.name, '>', valid_dimensions, '>', self.name
            arg.addClient(self)
        valid_dimensions.sort()
        self.valid_dimensions = tuple(valid_dimensions)
        self.args = args
        if 'name' in kw:
            self.name = kw.pop('name')
        self.setup(**kw)
    
    def setup(self):
        pass
    
    def update(self):
        for scope_t in self.assignments:
            scope = dict(zip(self.valid_dimensions, scope_t))
            input_nums = [arg.outputOrdinalFor(scope) for arg in self.args]
            self.assignments[scope_t] = tuple(input_nums)
        self._update_from_assignments()
        calc = self.makeCalcFunction()
        self.values = [nullor(self.name, calc, self.recycling)(*[a.values[i] for (i,a) in zip(u, self.args)]) for u in self.uniq]
    

class _LeafDefn(_Defn):
    """An input to the calculator, ie: a Defn with no inputs itself.
    
    This class is incomplete - subclasses provide:
        makeDefaultSetting()
        adaptSetting(setting)
        makeCells(input_soup)"""
    
    args = ()
    name = None
    name_required = True
    
    # This can be overriden in a subclass or the constuctor.
    valid_dimensions = ()
    
    array_template = None
    internal_dimensions = ()
    
    def __init__(self, name=None, extra_label=None,
            dimensions=None, independent_by_default=None):
        _Defn.__init__(self)
        if dimensions is not None:
            assert type(dimensions) in [list, tuple], type(dimensions)
            self.valid_dimensions = tuple(dimensions)
        if independent_by_default is not None:
            self.independent_by_default = independent_by_default
        if name is not None:
            self.name = name
        if self.name_required:
            assert isinstance(self.name, basestring), self.name
        if extra_label is not None:
            self.name = self.name + extra_label
    
    def getDefaultSetting(self):
        if (getattr(self, '_default_setting', None) is None or
                self.independent_by_default):
            self._default_setting = self.makeDefaultSetting()
        return self._default_setting
    
    def update(self):
        self._update_from_assignments()
        gdv = lambda x:x.getDefaultValue()
        self.values = [nullor(self.name, gdv)(u) for u in self.uniq]
    
    def assign(self, settings):
        for (scope, setting) in settings:
            self.checkSettingIsValid(setting)
            for scope_t in scope:
                assert scope_t in self.assignments, scope_t
                self.assignments[scope_t] = setting
    
    #def getNumFreeParams(self):
    #    return sum(setting.getNumFreeParams() for setting in self.uniq)
    
    def getAllDefaultValues(self, scope):
        return [self.assignments[s].getDefaultValue() for s in scope]
    
    def getCurrentBounds(self, scope):
        lowest = highest = None
        for s in scope:
            (lower, init, upper)  = self.assignments[s].getBounds()
            if upper == lower: continue
            if lowest is None or lower < lowest: lowest = lower
            if highest is None or upper > highest: highest = upper
        
        if lowest is None or highest is None:
            # All current settings are consts so use the class defaults
            (lowest, default, highest) = self.getDefaultSetting().getBounds()
        return (lowest, highest)
    
    def __repr__(self):
        return "%s(%s)" % (self.__class__.__name__,
                self._local_repr(col_width=6, max_width=60))
    
    def _local_repr(self, col_width, max_width):
        template = "%%%s.%sf" % (col_width, (col_width-1)//2)
        assignments = []
        for (i,a) in self.assignments.items():
            if a is None:
                assignments.append('None')
            elif a.is_const:
                if isinstance(a.value, float):
                    assignments.append(template % a.value)
                else:
                    assignments.append(a.value)
            else:
                assignments.append('Var') # %s' % str(i))
        return '%-20s%s' % (self.name[:19],
                _fmtrow(col_width+1, assignments, max_width))
    

class _ParameterController(object):
    """Holds a set of activated CalculationDefns, including their parameter
    scopes.  Makes calculators on demand."""
    
    def __init__(self, top_defn):
        # topological sort
        indegree = {id(top_defn):0}
        Q = [top_defn]
        while Q:
            pd = Q.pop(0)
            for arg in pd.args:
                arg_id = id(arg)
                if arg_id in indegree:
                    indegree[arg_id] += 1
                else:
                    indegree[arg_id] = 1
                    Q.append(arg)
        topdown = []
        Q = [top_defn]
        while Q:
            pd = Q.pop(0)
            topdown.append(pd)
            for arg in pd.args:
                arg_id = id(arg)
                indegree[arg_id] -= 1
                if indegree[arg_id] == 0:
                    Q.append(arg)
        
        # propagate categories downwards
        top_defn.assignments = {}
        for pd in topdown:
            pd.assignments = {}
            for client in pd.clients:
                scopes = client.getRequiredScopes(pd.valid_dimensions)
                # print pd.valid_dimensions, pd.name, '<', scopes, '<', client.name, client.valid_dimensions
                pd.addScopes(scopes)
            if not pd.assignments:
                pd.addScopes([{}])
            pd.activated = True
        
        self.defns = topdown[::-1]
        
        self.defn_for = {}
        for defn in self.defns:
            #if not defn.args:
            #assert defn.name not in self.defn_for, defn.name
            if defn.name in self.defn_for:
                self.defn_for[defn.name] = None
                # duplicate
            else:
                self.defn_for[defn.name] = defn
        
        self.update(self.defns)
    
    def getParamNames(self, scalar_only=False):
        return [defn.name for defn in self.defns if defn.user_param and
            (defn.numeric or not scalar_only)]
    
    def getUsedDimensions(self, par_name):
        return self.defn_for[par_name].usedDimensions()
    
    def getParamValue(self, par_name, *args, **kw):
        callback = self._makeValueCallback(None, None)
        defn = self.defn_for[par_name]
        posn = defn._getPosnForScope(*args, **kw)
        return callback(defn, posn)
    
    def getParamInterval(self, par_name, *args, **kw):
        dropoff = kw.pop('dropoff', None)
        p = kw.pop('p', None)
        if dropoff is None and p is None:
            p = 0.05
        callback = self._makeValueCallback(dropoff, p, kw.pop('xtol', None))
        defn = self.defn_for[par_name]
        posn = defn._getPosnForScope(*args, **kw)
        return callback(defn, posn)
    
    def getFinalResult(self):
        return self.defns[-1].getCurrentValueForScope()
    
    def getParamValueDict(self, dimensions, p=None, dropoff=None,
            params=None, xtol=None):
        callback = self._makeValueCallback(dropoff, p, xtol)
        if params is None:
            params = self.getParamNames(scalar_only=True)
        result = {}
        for param_name in params:
            ev = self.defn_for[param_name]
            ev.fillParValueDict(result, dimensions, callback)
        return result
    
    def _makeValueCallback(self, dropoff, p, xtol=None):
        """Make a setting -> value function"""
        if p is not None:
            assert dropoff is None, (p, dropoff)
            dropoff = chdtri(1, p) / 2.0
        if dropoff is None:
            def callback(defn, posn):
                return defn.values[posn]
        else:
            assert dropoff > 0, dropoff
            def callback(defn, posn):
                lc = self.makeCalculator(variable=defn.uniq[posn])
                assert len(lc.opt_pars) == 1, lc.opt_pars
                opt_par = lc.opt_pars[0]
                return lc._getCurrentCellInterval(opt_par, dropoff, xtol)
        return callback
    
    def update(self, changed=None):
        if changed is None:
            changed = self.defns # all
        # use topological sort order
        Q = set(id(defn) for defn in changed)
        for defn in self.defns:
            if id(defn) in Q:
                defn.update()
                for c in defn.clients:
                    Q.add(id(c))
    
    def __repr__(self):
        col_width = 6
        max_width = 60
        parts = [defn._local_repr(col_width, max_width) for defn in self.defns
                if not isinstance(defn, SelectFromDimension)]
        return '\n'.join(parts)
    
