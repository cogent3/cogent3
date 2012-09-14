#!/usr/bin/env python
from __future__ import division, with_statement
import warnings
import numpy
from contextlib import contextmanager
from .setting import Var, ConstVal
from .calculation import Calculator
from cogent.util import parallel
from cogent.maths.stats.distribution import chdtri
from cogent.maths.optimisers import MaximumEvaluationsReached

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.5.3"
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
        """Low level method for extracting values.  Pushes values of this 
        particular parameter/defn into the dict tree 'result', 
        eg: length_defn.fillParValueDict(['edge']) populates 'result' like
        {'length':{'human':1.0, 'mouse':1.0}}"""
        
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
    
    # These can be overriden in a subclass or the constuctor.
    valid_dimensions = ()
    numeric = False
    
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
    
    def assignAll(self, scope_spec=None, value=None,
            lower=None, upper=None, const=None, independent=None):
        settings = []
        if const is None:
            const = self.const_by_default
        
        for scope in self.interpretScopes(
                independent=independent, **(scope_spec or {})):
            if value is None:
                s_value = self.getMeanCurrentValue(scope)
            else:       
                s_value = self.unwrapValue(value)
            
            if const:
                setting = ConstVal(s_value)
            elif not self.numeric:
                if lower is not None or upper is not None:
                    raise ValueError(
                            "Non-scalar input '%s' doesn't support bounds"
                            % self.name)
                setting = Var((None, s_value, None))
            else:
                (s_lower, s_upper) = self.getCurrentBounds(scope)
                if lower is not None:
                    s_lower = lower
                if upper is not None:
                    s_upper = upper
                    
                if s_lower > s_upper:
                    raise ValueError("Bounds: upper < lower")
                elif (s_lower is not None) and s_value < s_lower:
                    s_value = s_lower
                    warnings.warn("Value of %s increased to keep within bounds" 
                                % self.name, stacklevel=3)
                elif (s_upper is not None) and s_value > s_upper:
                    s_value = s_upper
                    warnings.warn("Value of %s decreased to keep within bounds" 
                                % self.name, stacklevel=3)
                setting = Var((s_lower, s_value, s_upper))
            self.checkSettingIsValid(setting)
            settings.append((scope, setting))
            
        for (scope, setting) in settings:
            for scope_t in scope:
                assert scope_t in self.assignments, scope_t
                self.assignments[scope_t] = setting
    
    def getMeanCurrentValue(self, scope):
        values = [self.assignments[s].getDefaultValue() for s in scope]
        if len(values) == 1:
            s_value = values[0]
        else:
            s_value = sum(values) / len(values)
            for value in values:
                if not numpy.all(value==s_value):
                    warnings.warn("Used mean of %s %s values" % 
                            (len(values), self.name), stacklevel=4)
                    break
        return s_value
                    
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
            elif a.is_constant:
                if isinstance(a.value, float):
                    assignments.append(template % a.value)
                else:
                    assignments.append(a.value)
            else:
                assignments.append('Var') # %s' % str(i))
        return '%-20s%s' % (self.name[:19],
                _fmtrow(col_width+1, assignments, max_width))
    

class ParameterController(object):
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
        
        self._changed = set()
        self._update_suspended = False
        self.updateIntermediateValues(self.defns)
        self.setupParallelContext()
    
    def getParamNames(self, scalar_only=False):
        """The names of the numerical inputs to the calculation."""
        return [defn.name for defn in self.defns if defn.user_param and
            (defn.numeric or not scalar_only)]
    
    def getUsedDimensions(self, par_name):
        return self.defn_for[par_name].usedDimensions()
    
    def getParamValue(self, par_name, *args, **kw):
        """The value for 'par_name'.  Additional arguments specify the scope.
        Despite the name intermediate values can also be retrieved this way."""
        callback = self._makeValueCallback(None, None)
        defn = self.defn_for[par_name]
        posn = defn._getPosnForScope(*args, **kw)
        return callback(defn, posn)
    
    def getParamInterval(self, par_name, *args, **kw):
        """Confidence interval for 'par_name' found by adjusting the 
        single parameter until the final result falls by 'dropoff', which
        can be specified directly or via 'p' as chdtri(1, p).  Additional
        arguments are taken to specify the scope."""
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
        """A dict tree of parameter values, with parameter names as the
        top level keys, and the various dimensions ('edge', 'bin', etc.)
        supplying lower level keys: edge names, bin names etc.
        If 'p' or 'dropoff' is specified returns chi-square intervals instead
        of simple values."""
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
    
    @contextmanager
    def updatesPostponed(self):
        "Temporarily turn off calculation for faster input setting"
        (old, self._update_suspended) = (self._update_suspended, True)
        yield
        self._update_suspended = old
        self._updateIntermediateValues()
        
    def updateIntermediateValues(self, changed=None):
        if changed is None:
            changed = self.defns # all
        self._changed.update(id(defn) for defn in changed)
        self._updateIntermediateValues()
    
    def _updateIntermediateValues(self):
        if self._update_suspended:
            return
        # use topological sort order
        # xxx parallel context check?
        for defn in self.defns:
            if id(defn) in self._changed:
                defn.update()
                for c in defn.clients:
                    self._changed.add(id(c))
        self._changed.clear()
    
    def assignAll(self, par_name, *args, **kw):
        defn = self.defn_for[par_name]
        if not isinstance(defn, _LeafDefn):
            args = ' and '.join(['"%s"' % a.name for a in defn.args])
            msg = '"%s" is not settable as it is derived from %s.' % (
                    par_name, args)
            raise ValueError(msg)
        defn.assignAll(*args, **kw)
        self.updateIntermediateValues([defn])
        
    def measureEvalsPerSecond(self, *args, **kw):
        return self.makeCalculator().measureEvalsPerSecond(*args, **kw)
    
    def setupParallelContext(self, parallel_split=None):
        self.overall_parallel_context = parallel.getContext()
        with parallel.split(parallel_split) as parallel_context:
            parallel_context = parallel_context.getCommunicator()
            self.remaining_parallel_context = parallel.getContext()
            if 'parallel_context' in self.defn_for:
                self.assignAll(
                    'parallel_context', value=parallel_context, const=True)
    
    def makeCalculator(self, calculatorClass=None, variable=None, **kw):
        cells = []
        input_soup = {}
        for defn in self.defns:
            defn.update()
            (newcells, outputs) = defn.makeCells(input_soup, variable)
            cells.extend(newcells)
            input_soup[id(defn)] = outputs
        if calculatorClass is None:
            calculatorClass = Calculator
        kw['overall_parallel_context'] = self.overall_parallel_context
        kw['remaining_parallel_context'] = self.remaining_parallel_context
        return calculatorClass(cells, input_soup, **kw)
    
    def updateFromCalculator(self, calc):
        changed = []
        for defn in self.defn_for.values():
            if isinstance(defn, _LeafDefn):
                defn.updateFromCalculator(calc)
                changed.append(defn)
        self.updateIntermediateValues(changed)
    
    def getNumFreeParams(self):
        return sum(defn.getNumFreeParams() for defn in self.defns if isinstance(defn, _LeafDefn))

    def optimise(self, local=None, 
            filename=None, interval=None,
            limit_action='warn',  max_evaluations=None, 
            tolerance=1e-6, global_tolerance=1e-1, **kw):
        """Find input values that optimise this function.
        'local' controls the choice of optimiser, the default being to run
        both the global and local optimisers. 'filename' and 'interval'
        control checkpointing.  Unknown keyword arguments get passed on to
        the optimiser(s)."""
        return_calculator = kw.pop('return_calculator', False) # only for debug
        for n in ['local', 'filename', 'interval', 'max_evaluations', 
                'tolerance', 'global_tolerance']:
            kw[n] = locals()[n]
        lc = self.makeCalculator()
        try:
            lc.optimise(**kw)
        except MaximumEvaluationsReached, detail:
            evals = detail[0]
            err_msg = 'FORCED EXIT from optimiser after %s evaluations' % evals
            if limit_action == 'ignore':
                pass
            elif limit_action == 'warn':
                warnings.warn(err_msg, stacklevel=2)
            else:
                raise ArithmeticError(err_msg)    
        finally:
            self.updateFromCalculator(lc)
        if return_calculator:
            return lc
    
    def graphviz(self, **kw):
        lc = self.makeCalculator()
        return lc.graphviz(**kw)
        


