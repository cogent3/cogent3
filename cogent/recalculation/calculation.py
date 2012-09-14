#!/usr/bin/env python
from __future__ import division, with_statement
import numpy
Float = numpy.core.numerictypes.sctype2char(float)
import time, warnings
from cogent.maths.solve import find_root
from cogent.util import parallel
from cogent.maths.optimisers import maximise, ParameterOutOfBoundsError


import os
TRACE_DEFAULT = os.environ.has_key('COGENT_TRACE')
TRACE_SCALE = 100000

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Peter Maxwell", "Gavin Huttley", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"

# This is the 'live' layer of the recalculation system
# Cells and OptPars are held by a Calculator
# For docstring see definitions.py

class CalculationInterupted(Exception):
    pass

class OptPar(object):
    """One parameter, as seen by the optimiser, eg: length of one edge.
    An OptPar reports changes to the ParameterValueSet for its parameter.
    """
    
    is_constant = False
    recycled = False
    args = ()
    # Use of __slots__ here and in Cell gives 8% speedup on small calculators.
    __slots__ = ['clients', 'client_ranks', 'name', 'lower', 'default_value',
            'upper', 'scope', 'order', 'label', 'consequences', 'rank']
    
    def __init__(self, name, scope, bounds):
        self.clients = []
        self.client_ranks = []
        self.name = name
        for (attr, v) in zip(['lower', 'default_value', 'upper'], bounds):
            setattr(self, attr, float(v))
        
        # controls order in optimiser - group for LF
        self.scope = scope
        self.order = (len(scope), scope and min(scope), name)
        self.label = self.name
    
    def addClient(self, client):
        self.clients.append(client)
    
    def __cmp__(self, other):
        # optimisation is more efficient if params for one edge are neighbours
        return cmp(self.order, other.order)
    
    def __repr__(self):
        return '%s(%s)' % (self.__class__.__name__, self.label)
    
    def getOptimiserBounds(self):
        lower = self.transformToOptimiser(self.lower)
        upper = self.transformToOptimiser(self.upper)
        return (lower, upper)
    
    def transformFromOptimiser(self, value):
        return value
    
    def transformToOptimiser(self, value):
        return value
    

class LogOptPar(OptPar):
    # For ratios, optimiser sees log(param value).  Conversions to/from
    # optimiser representation are only done by Calculator.change(),
    # .getValueArray() and .getBoundsArrrays().
    
    def transformFromOptimiser(self, value):
        return numpy.exp(value)
    
    def transformToOptimiser(self, value):
        try:
            return numpy.log(value)
        except OverflowError:
            raise OverflowError('log(%s)' % value)
    

class EvaluatedCell(object):
    __slots__ = ['client_ranks', 'rank', 'calc', 'args', 'is_constant',
        'clients', 'failure_count', 'name', 'arg_ranks',
        'consequences', 'recycled', 'default']
    
    def __init__(self, name, calc, args, recycling=None, default=None):
        self.name = name
        self.rank = None
        self.calc = calc
        self.default = default
        self.args = tuple(args)
        
        self.recycled = recycling
        if recycling:
            self.args = (self,) + self.args
        
        self.is_constant = True
        for arg in args:
            arg.addClient(self)
            if not arg.is_constant:
                self.is_constant = False
        
        self.clients = []
        self.client_ranks = []
        self.failure_count = 0
    
    def addClient(self, client):
        self.clients.append(client)
    
    def update(self, data):
        data[self.rank] = self.calc(
                *[data[arg_rank] for arg_rank in self.arg_ranks])
    
    def prime(self, data_sets):
        if self.is_constant:
            # Just calc once
            self.update(data_sets[0])
            for data in data_sets[1:]:
                data[self.rank] = data_sets[0][self.rank]
        else:
            for data in data_sets:
                self.update(data)
    
    def reportError(self, detail, data):
        self.failure_count += 1
        if self.failure_count <= 5:
            print ("%s in calculating %s:",
                    detail.__class__.__name__, self.name)
        if self.failure_count == 5:
            print "Additional failures of this type will not be reported."
        if self.failure_count < 2:
            print '%s inputs were:', len(self.arg_ranks)
            for (i, arg) in enumerate(self.arg_ranks):
                print '%s: ' % i + repr(data[arg])
    

class ConstCell(object):
    __slots__ = ['name', 'scope', 'value', 'rank', 'consequences', 'clients']
    
    recycled = False
    is_constant = True
    args = ()
    
    def __init__(self, name, value):
        self.name = name
        self.clients = []
        self.value = value
    
    def addClient(self, client):
        self.clients.append(client)
    

class Calculator(object):
    """A complete hierarchical function with N evaluation steps to call
    for each change of inputs.  Made by a ParameterController."""
    
    def __init__(self, cells, defns, remaining_parallel_context=None,
                overall_parallel_context=None, trace=None, with_undo=True):
        if trace is None:
            trace = TRACE_DEFAULT
        self.overall_parallel_context = overall_parallel_context
        self.remaining_parallel_context = remaining_parallel_context
        self.with_undo = with_undo
        self.results_by_id = defns
        self.opt_pars = []
        other_cells = []
        for cell in cells:
            if isinstance(cell, OptPar):
                self.opt_pars.append(cell)
            else:
                other_cells.append(cell)
        self._cells = self.opt_pars + other_cells
        data_sets = [[0], [0,1]][self.with_undo]
        self.cell_values = [[None]*len(self._cells) for switch in data_sets]
        self.arg_ranks = [[] for cell in self._cells]
        for (i, cell) in enumerate(self._cells):
            cell.rank = i
            cell.consequences = {}
            if isinstance(cell, OptPar):
                for switch in data_sets:
                    self.cell_values[switch][i] = cell.default_value
            elif isinstance(cell, ConstCell):
                for switch in data_sets:
                    self.cell_values[switch][i] = cell.value
            elif isinstance(cell, EvaluatedCell):
                cell.arg_ranks = []
                for arg in cell.args:
                    if hasattr(arg, 'client_ranks'):
                        arg.client_ranks.append(i)
                    self.arg_ranks[i].append(arg.rank)
                    cell.arg_ranks.append(arg.rank)
                
                with parallel.parallel_context(self.remaining_parallel_context):
                    try:
                        cell.prime(self.cell_values)
                    except KeyboardInterrupt:
                        raise
                    except Exception, detail:
                        print ("Failed initial calculation of %s"
                                % cell.name)
                        raise
            else:
                raise RuntimeError('Unexpected Cell type %s' % type(cell))
        
        self._switch = 0
        self.recycled_cells = [
                cell.rank for cell in self._cells if cell.recycled]
        self.spare = [None] * len (self._cells)
        
        for cell in self._cells[::-1]:
            for arg in cell.args:
                arg.consequences[cell.rank] = True
                arg.consequences.update(cell.consequences)
        
        self._programs = {}
        # Just for timings pre-calc these
        for opt_par in self.opt_pars:
            self.cellsChangedBy([(opt_par.rank, None)])
        
        self.last_values = self.getValueArray()
        self.last_undo = []
        self.elapsed_time = 0.0
        self.evaluations = 0
        self.setTracing(trace)
        self.optimised = False
    
    def _graphviz(self):
        """A string in the 'dot' graph description language used by the
        program 'Graphviz'.  One box per cell, grouped by Defn."""
        
        lines = ['digraph G {\n rankdir = LR\n ranksep = 1\n']
        evs = []
        for cell in self._cells:
            if cell.name not in evs:
                evs.append(cell.name)
        nodes = dict([(name, []) for name in evs])
        edges = []
        for cell in self._cells:
            if hasattr(cell, 'name'):
                nodes[cell.name].append(cell)
                for arg in cell.args:
                    if arg is not cell:
                        edges.append('"%s":%s -> "%s":%s' %
                                (arg.name, arg.rank, cell.name, cell.rank))
        for name in evs:
            all_const = True
            some_const = False
            enodes = [name.replace('edge', 'QQQ')]
            for cell in nodes[name]:
                value = self._getCurrentCellValue(cell)
                if isinstance(value, float):
                    label = '%5.2e' % value
                else:
                    label = '[]'
                label = '<%s> %s' % (cell.rank, label)
                enodes.append(label)
                all_const = all_const and cell.is_constant
                some_const = some_const or cell.is_constant
            enodes = '|'.join(enodes)
            colour = ['', ' fillcolor=gray90, style=filled,'][some_const]
            colour = [colour, ' fillcolor=gray, style=filled,'][all_const]
            lines.append('"%s" [shape = "record",%s label="%s"];' %
                    (name, colour, enodes))
        lines.extend(edges)
        lines.append('}')
        return '\n'.join(lines).replace('edge', 'egde').replace('QQQ', 'edge')
    
    def graphviz(self, keep=False):
        """Use Graphviz to display a graph representing the inner workings of
        the calculator.  Leaves behind a temporary file (so that Graphviz can
        redraw it with different settings) unless 'keep' is False"""
        
        import tempfile, os, sys
        
        if sys.platform != 'darwin':
            raise NotImplementedError, "Graphviz support Mac only at present"
        
        GRAPHVIZ = '/Applications/Graphviz.app'
        # test that graphviz is installed
        if not os.path.exists(GRAPHVIZ):
            raise RuntimeError('%s not present' % GRAPHVIZ)
        
        text = self._graphviz()
        
        fn = tempfile.mktemp(prefix="calc_", suffix=".dot")
        f = open(fn, 'w')
        f.write(text)
        f.close()
        # Mac specific!
        # Specify Graphviz as ".dot" can mean other things.
        # Would be sensible to eventually use LaunchServices.
        os.system('open -a "%s" "%s"' % (GRAPHVIZ, fn))
        if not keep:
            time.sleep(5)
            os.remove(fn)
    
    def optimise(self, **kw):
        x = self.getValueArray()
        bounds = self.getBoundsVectors()
        maximise(self, x, bounds, **kw)
        self.optimised = True
    
    def setTracing(self, trace=False):
        """With 'trace' true every evaluated is printed.  Useful for profiling
        and debugging."""
        
        self.trace = trace
        if trace:
            print
            n_opars = len(self.opt_pars)
            n_cells = len([c for c in self._cells if not c.is_constant])
            print n_opars, "OptPars and", n_cells - n_opars, "derived values"
            print 'OptPars: ', ', '.join([par.name for par in self.opt_pars])
            print "Times in 1/%sths of a second" % TRACE_SCALE

            groups = []
            groupd = {}
            for cell in self._cells:
                if cell.is_constant or not isinstance(cell, EvaluatedCell):
                    continue
                if cell.name not in groupd:
                    group = []
                    groups.append((cell.name, group))
                    groupd[cell.name] = group
                groupd[cell.name].append(cell)
            
            widths = []
            for (name, cells) in groups:
                width = 4 + len(cells)
                widths.append(min(15, width))
            self._cellsGroupedForDisplay = zip(groups, widths)
            for ((name, cells), width) in self._cellsGroupedForDisplay:
                print name[:width].ljust(width), '|',
            print
            for width in widths:
                print '-' * width, '|',
            print
    
    def getValueArray(self):
        """This being a caching function, you can ask it for its current
        input!  Handy for initialising the optimiser."""
        values = [p.transformToOptimiser(self._getCurrentCellValue(p))
                for p in self.opt_pars]
        return values
    
    # getBoundsVectors and testoptparvector make up the old LikelihoodFunction
    # interface expected by the optimiser.
    
    def getBoundsVectors(self):
        """2 arrays: minimums, maximums"""
        lower = numpy.zeros([len(self.opt_pars)], Float)
        upper = numpy.zeros([len(self.opt_pars)], Float)
        for (i, opt_par) in enumerate(self.opt_pars):
            (lb, ub) = opt_par.getOptimiserBounds()
            lower[i] = lb
            upper[i] = ub
        return (lower, upper)
    
    def fuzz(self, random_series=None, seed=None):
        # Slight randomisation suitable for removing right-on-the-
        # ridge starting points before local optimisation.
        if random_series is None:
            import random
            random_series = random.Random()
        if seed is not None:
            random_series.seed(seed)
        X = self.getValueArray()
        for (i, (l,u)) in enumerate(zip(*self.getBoundsVectors())):
            sign = random_series.choice([-1, +1])
            step = random_series.uniform(+0.05, +0.025)
            X[i] = max(l,min(u,(1.0 + sign*step*X[i])))
        self.testoptparvector(X)
        self.optimised = False
    
    def testoptparvector(self, values):
        """AKA self().  Called by optimisers.  Returns the output value
        after doing any recalculation required for the new input 'values'
        array"""
        
        assert len(values) == len(self.opt_pars)
        changes = [(i, new) for (i, (old, new))
                in enumerate(zip(self.last_values, values))
                if old != new]
        return self.change(changes)
    
    __call__ = testoptparvector
    
    def testfunction(self):
        """Return the current output value without changing any inputs"""
        return self._getCurrentCellValue(self._cells[-1])
    
    def change(self, changes):
        """Returns the output value after applying 'changes', a list of
        (optimisable_parameter_ordinal, new_value) tuples."""
        
        t0 = time.time()
        self.evaluations += 1
        assert parallel.getContext() is self.overall_parallel_context, (
            parallel.getContext(), self.overall_parallel_context)
        
        # If ALL of the changes made in the last step are reversed in this step
        # then it is safe to undo them first, taking advantage of the 1-deep
        # cache.
        if self.with_undo and self.last_undo:
            for (i, v) in self.last_undo:
                if (i,v) not in changes:
                    break
            else:
                changes = [ch for ch in changes if ch not in self.last_undo]
                self._switch = not self._switch
                for (i, v) in self.last_undo:
                    self.last_values[i] = v
        
        self.last_undo = []
        program = self.cellsChangedBy(changes)
        
        if self.with_undo:
            self._switch = not self._switch
            data = self.cell_values[self._switch]
            base = self.cell_values[not self._switch]
            
            # recycle and undo interact in bad ways
            for rank in self.recycled_cells:
                if data[rank] is not base[rank]:
                    self.spare[rank] = data[rank]
            data[:] = base[:]
            for cell in program:
                if cell.recycled:
                    if data[cell.rank] is base[cell.rank]:
                        data[cell.rank]=self.spare[cell.rank]
                        assert data[cell.rank] is not base[cell.rank]
        else:
            data = self.cell_values[self._switch]
        
        # Set new OptPar values
        changed_optpars = []
        for (i, v)  in changes:
            if i < len(self.opt_pars):
                assert isinstance(v*1.0, float), v
                changed_optpars.append((i, self.last_values[i]))
                self.last_values[i] = v
                data[i] = self.opt_pars[i].transformFromOptimiser(v)
            else:
                data[i] = v
        
        with parallel.parallel_context(self.remaining_parallel_context):
            try:
                if self.trace:
                    self.tracingUpdate(changes, program, data)
                else:
                    self.plainUpdate(program, data)
                
                # if non-optimiser parameter was set then undo is invalid
                if (self.last_undo and
                        max(self.last_undo)[0] >= len(self.opt_pars)):
                    self.last_undo = []
                else:
                    self.last_undo = changed_optpars
            
            except CalculationInterupted, detail:
                if self.with_undo:
                    self._switch = not self._switch
                for (i,v) in changed_optpars:
                    self.last_values[i] = v
                self.last_undo = []
                (cell, exception) = detail.args
                raise exception
            
            finally:
                self.elapsed_time += time.time() - t0
        
        return self.cell_values[self._switch][-1]
    
    def cellsChangedBy(self, changes):
        # What OptPars have been changed determines cells to update
        change_key = dict(changes).keys()
        change_key.sort()
        change_key = tuple(change_key)
        if change_key in self._programs:
            program = self._programs[change_key]
        else:
            # Make a list of the cells to update and cache it.
            consequences = {}
            for i in change_key:
                consequences.update(self._cells[i].consequences)
            self._programs[change_key] = program = (
                [cell for cell in self._cells if cell.rank in consequences])
        return program
    
    def plainUpdate(self, program, data):
        try:
            for cell in program:
                data[cell.rank] = cell.calc(*[data[a] for a in cell.arg_ranks])
        except ParameterOutOfBoundsError, detail:
            # Non-fatal error, just cancel this calculation.
            raise CalculationInterupted(cell, detail)
        except ArithmeticError, detail:
            # Non-fatal but unexpected error. Warn and cancel this calculation.
            cell.reportError(detail, data)
            raise CalculationInterupted(cell, detail)
    
    def tracingUpdate(self, changes, program, data):
        # Does the same thing as plainUpdate, but also produces lots of
        # output showing how long each step of the calculation takes.
        # One line per call, '-' for undo, '+' for calculation
                
        exception = None
        elapsed = {}
        for cell in program:
            try:
                t0 = time.time()
                data[cell.rank] = cell.calc(*[data[a] for a in cell.arg_ranks])
                t1 = time.time()
            except (ParameterOutOfBoundsError, ArithmeticError), exception:
                error_cell = cell
                break
            elapsed[cell.rank] = (t1-t0)
        
        tds = []
        for ((name, cells), width) in self._cellsGroupedForDisplay:
            text = ''.join([' +'[cell.rank in elapsed] for cell in cells])
            elap = sum([elapsed.get(cell.rank, 0) for cell in cells])
            if len(text) > width-4:
                edge_width = min(len(text), (width - 4 - 3)) // 2
                elipsis = ['   ','...'][not not text.strip()]
                text = text[:edge_width] + elipsis + text[-edge_width:]
            tds.append('%s%4s' % (text, int(TRACE_SCALE*elap+0.5) or ''))
        
        par_descs = []
        for (i,v) in changes:
            cell = self._cells[i]
            if isinstance(cell, OptPar):
                par_descs.append('%s=%8.6f' % (cell.name, v))
            else:
                par_descs.append('%s=?' % cell.name)
        par_descs = ', '.join(par_descs)[:22].ljust(22)
        print ' | '.join(tds+['']),
        if exception:
            print '%15s | %s' % ('', par_descs)
            error_cell.reportError(exception, data)
            raise CalculationInterupted(cell, exception)
        else:
            print '%-15s | %s' % (repr(data[-1])[:15], par_descs)
    
    def measureEvalsPerSecond(self, time_limit=1.0, wall=True, sa=False):
        # Returns an estimate of the number of evaluations per second
        # an each-optpar-in-turn simulated annealing type optimiser
        # can achive, spending not much more than 'time_limit' doing
        # so.  'wall'=False causes process time to be used instead of
        # wall time.
        # 'sa' makes it simulated-annealing-like, with frequent backtracks
        if wall:
            now = time.time
        else:
            now = time.clock
        x = self.getValueArray()
        samples = []
        elapsed = 0.0
        rounds_per_sample = 2
        comm = parallel.getCommunicator()
        while elapsed < time_limit and len(samples) < 5:
            time.sleep(0.01)
            t0 = now()
            last = []
            for j in range(rounds_per_sample):
                for (i,v) in enumerate(x):
                     # Not a real change, but works like one.
                    self.change(last + [(i, v)])
                    if sa and (i+j) % 2:
                        last = [(i, v)]
                    else:
                        last = []
            # Use one agreed on delta otherwise different cpus will finish the
            # loop at different times causing chaos.
            delta = comm.allreduce(now()-t0, parallel.MPI.MAX)
            if delta < 0.1:
                # time.clock is low res, so need to ensure each sample
                # is long enough to take SOME time.
                rounds_per_sample *= 2
                continue
            else:
                rate = rounds_per_sample * len(x) / delta
                samples.append(rate)
                elapsed += delta
        
        if wall:
            samples.sort()
            return samples[len(samples)//2]
        else:
            return sum(samples) / len(samples)
    
    def _getCurrentCellValue(self, cell):
        return self.cell_values[self._switch][cell.rank]
    
    def getCurrentCellValuesForDefn(self, defn):
        cells = self.results_by_id[id(defn)]
        return [self.cell_values[self._switch][cell.rank] for cell in cells]

    def __getBoundedRoot(self, func, origX, direction, bound, xtol):
        return find_root(func, origX, direction, bound, xtol=xtol,
                expected_exception = (
                    ParameterOutOfBoundsError, ArithmeticError))
    
    def _getCurrentCellInterval(self, opt_par, dropoff, xtol=None):
        # (min, opt, max) tuples for each parameter where f(min) ==
        # f(max) == f(opt)-dropoff.  Uses None when a bound is hit.
        #assert self.optimised, "Call optimise() first"
        origY = self.testfunction()
        (lower, upper) = opt_par.getOptimiserBounds()
        opt_value = self._getCurrentCellValue(opt_par)
        origX = opt_par.transformToOptimiser(opt_value)
        def func(x):
            Y = self.change([(opt_par.rank, x)])
            return Y - (origY - dropoff)
        try:
            lowX = self.__getBoundedRoot(func, origX, -1, lower, xtol)
            highX = self.__getBoundedRoot(func, origX, +1, upper, xtol)
        finally:
            func(origX)
        
        triple = []
        for x in [lowX, origX, highX]:
            if x is not None:
                x = opt_par.transformFromOptimiser(x)
            triple.append(x)
        return tuple(triple)
        
    
