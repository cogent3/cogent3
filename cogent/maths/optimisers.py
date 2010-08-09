#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Local or Global-then-local optimisation with progress display
"""

from cogent.util import progress_display as UI
from simannealingoptimiser import SimulatedAnnealing
from scipy_optimisers import DownhillSimplex, Powell
import warnings

GlobalOptimiser = SimulatedAnnealing
LocalOptimiser = Powell

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Andrew Butterfield", "Peter Maxwell", "Gavin Huttley"]
__license__ = "GPL"
__version__ = "1.5.0.dev"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"

def unsteadyProgressIndicator(display_progress, label='', start=0.0, end=1.0):
    template = 'f = % #10.6g  ±  % 9.3e   evals = %6i '
    label = label.rjust(5)
    goal = [1.0e-20]
    def _display_progress(remaining, *args):
        if remaining > goal[0]:
            goal[0] = remaining
        progress = (goal[0]-remaining)/goal[0] * (end-start) + start
        msg = template % args + label
        return display_progress(msg, progress=progress, current=0)
    return _display_progress

@UI.display_wrap
def optimise(f, x, bounds, local=None, filename=None, interval=None,
        max_restarts=None, max_evaluations=None,
        tolerance=1e-6, global_tolerance=1e-1, ui=None, **kw):
    """Find input values that optimise this function.
    'local' controls the choice of optimiser, the default being to run
    both the global and local optimisers. 'filename' and 'interval'
    control checkpointing.  Unknown keyword arguments get passed on to
    the optimiser(s)."""
    do_global = (not local) or local is None
    do_local = local or local is None
        
    # Global optimisation
    if do_global:
        if not do_local:
            wanings.warn(
                'local=False causes the post-global optimisation local '
                '"polishing" optimisation to be skipped entirely, which seems '
                'pointless, so its meaning may change to a simple boolean '
                'flag: local or global.')
            # It also needlessly complicates this function.
            gend = 1.0
        else:
            gend = 0.9
        callback = unsteadyProgressIndicator(ui.display, 'Global', 0.0, gend)
        gtol = [tolerance, global_tolerance][do_local]
        opt = GlobalOptimiser(f, x, bounds, tolerance=gtol,
                max_evaluations=max_evaluations, **kw)
        opt.setCheckpointing(filename=filename, interval=interval)
        result = opt.run(callback)
    else:
        gend = 0.0
        for k in kw:
            warnings.warn('Unused arg for local alignment: ' + k)
    
    # Local optimisation
    if do_local:
        callback = unsteadyProgressIndicator(ui.display, 'Local', gend, 1.0)
        #ui.display('local opt', 1.0-per_opt, per_opt)
        opt = LocalOptimiser(f, x, bounds, tolerance=tolerance,
            max_restarts=max_restarts, max_evaluations=max_evaluations)
        result = opt.run(callback)
        
    return result

