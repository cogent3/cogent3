#!/usr/bin/env python
from distutils.core import setup, Command
from distutils.extension import Extension
import sys, os, re, subprocess

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2011, The Cogent Project"
__contributors__ = ["Peter Maxwell", "Gavin Huttley", "Matthew Wakefield",
                    "Greg Caporaso", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"

# Check Python version, no point installing if unsupported version inplace
if sys.version_info < (2, 6):
    py_version = ".".join([str(n) for n in sys.version_info])
    raise RuntimeError("Python-2.6 or greater is required, Python-%s used." % py_version)


# Check Numpy version, no point installing if unsupported version inplace
try:
    import numpy
except ImportError:
    raise RuntimeError("Numpy required but not found.")

numpy_version = re.split("[^\d]", numpy.__version__)
numpy_version_info = tuple([int(i) for i in numpy_version if i.isdigit()])
if numpy_version_info < (1, 3):
    raise RuntimeError("Numpy-1.3 is required, %s found." % numpy_version)

# Find arrayobject.h on any system
numpy_include_dir = numpy.get_include()


# On windows with no commandline probably means we want to build an installer.
if sys.platform == "win32" and len(sys.argv) < 2:
    sys.argv[1:] = ["bdist_wininst"]


# A new command for predist, ie: pyrexc but no compile.
class NullCommand(Command):
    description = "Generate .c files from .pyx files"
    # List of option tuples: long name, short name (or None), and help string.
    user_options = [] #[('', '', ""),]
    def initialize_options (self):
        pass
    def finalize_options (self):
        pass
    def run (self):
        pass
        
class BuildDocumentation(NullCommand):
    description = "Generate HTML documentation and .c files"
    def run (self):
        # Restructured Text -> HTML
        try:
            import sphinx
        except ImportError:
            print "Failed to build html due to ImportErrors for sphinx"
            return
        cwd = os.getcwd()
        os.chdir('doc')
        subprocess.call(["make", "html"])
        os.chdir(cwd)
        print "Built index.html"

# Cython is now run via the Cythonize function rather than monkeypatched into 
# distutils, so these legacy commands don't need to do anything extra.
extra_commands = {
    'pyrexc': NullCommand,
    'cython': NullCommand,
    'predist': BuildDocumentation}


# Compiling Pyrex modules to .c and .so, if possible and necessary
try:
    if 'DONT_USE_PYREX' in os.environ:
        raise ImportError
    from Cython.Compiler.Version import version
    version = tuple([int(v) \
        for v in re.split("[^\d]", version) if v.isdigit()])
    if version < (0, 17, 1):
        print "Your Cython version is too old"
        raise ImportError
except ImportError:
    source_suffix = '.c'
    cythonize = lambda x:x
    print "No Cython, will compile from .c files"
    for cmd in extra_commands:
        if cmd in sys.argv:
            print "'%s' command not available without Cython" % cmd
            sys.exit(1)
else:
    from Cython.Build import cythonize
    source_suffix = '.pyx'


# Save some repetitive typing.  We have all compiled modules in place
# with their python siblings, so their paths and names are the same.
def CythonExtension(module_name, **kw):
    path = module_name.replace('.', '/')
    return Extension(module_name, [path + source_suffix], **kw)


short_description = "COmparative GENomics Toolkit"

# This ends up displayed by the installer
long_description = """Cogent
A toolkit for statistical analysis of biological sequences.
Version %s.
""" % __version__

setup(
    name="cogent",
    version=__version__,
    url="http://github.com/pycogent/pycogent",
    author="Gavin Huttley, Rob Knight",
    author_email="gavin.huttley@anu.edu.au, rob@spot.colorado.edu",
    description=short_description,
    long_description=long_description,
    platforms=["any"],
    license=["GPL"],
    keywords=["biology", "genomics", "statistics", "phylogeny", "evolution",
                "bioinformatics"],
    classifiers=[
            "Development Status :: 5 - Production/Stable",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: GNU General Public License (GPL)",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            "Topic :: Software Development :: Libraries :: Python Modules",
            "Operating System :: OS Independent",
            ],
    packages=['cogent', 'cogent.align', 'cogent.align.weights', 'cogent.app',
                'cogent.cluster', 'cogent.core', 'cogent.data', 'cogent.db',
                'cogent.db.ensembl', 'cogent.draw',
                'cogent.evolve', 'cogent.format', 'cogent.maths',
                'cogent.maths.matrix', 'cogent.maths.stats',
                'cogent.maths.stats.cai', 'cogent.maths.unifrac',
                'cogent.maths.spatial', 'cogent.motif', 'cogent.parse',
                'cogent.phylo', 'cogent.recalculation', 'cogent.seqsim',
                'cogent.struct', 'cogent.util'],
    ext_modules=cythonize([
        CythonExtension("cogent.align._compare"),
        CythonExtension("cogent.align._pairwise_seqs"),
        CythonExtension("cogent.align._pairwise_pogs"),
        CythonExtension("cogent.evolve._solved_models"),
        CythonExtension("cogent.evolve._likelihood_tree"),
        CythonExtension("cogent.evolve._pairwise_distance"),
        CythonExtension("cogent.struct._asa"),
        CythonExtension("cogent.struct._contact"),
        CythonExtension("cogent.maths._period"),
        CythonExtension("cogent.maths.spatial.ckd3"),
    ]),
    include_dirs = [numpy_include_dir],
    cmdclass = extra_commands,
)
