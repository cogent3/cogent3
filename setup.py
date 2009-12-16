#!/usr/bin/env python
from distutils.core import setup
import sys, os, re, subprocess

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__contributors__ = ["Peter Maxwell", "Gavin Huttley", "Matthew Wakefield",
                    "Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"

doc_imports_failed = False
try:
    import sphinx
except ImportError:
    doc_imports_failed = True

# A new command for predist, ie: pyrexc but no compile.
import distutils.ccompiler
class NullCompiler(distutils.ccompiler.CCompiler):
    # this basically to stop pyrexc building binaries, just the .c files
    executables = ()
    def __init__(self, *args, **kw):
        pass
    def compile(self, *args, **kw):
        return []
    def link(self, *args, **kw):
        pass

# Pyrex makes some messy C code so limit some warnings when we know how.
import distutils.sysconfig
if (distutils.sysconfig.get_config_var('CC') or '').startswith("gcc"):
    pyrex_compile_options = ['-w']
else:
    pyrex_compile_options = []

# On windows with no commandline probably means we want to build an installer.
if sys.platform == "win32" and len(sys.argv) < 2:
    sys.argv[1:] = ["bdist_wininst"]

# Restructured Text -> HTML
def build_html():
    if doc_imports_failed:
        print "Failed to build html due to ImportErrors for sphinx"
        return
    cwd = os.getcwd()
    os.chdir('doc')
    subprocess.call(["make", "html"])
    os.chdir(cwd)
    print "Built index.html"

# Compiling Pyrex modules to .c and .so
include_path = os.path.join(os.getcwd(), 'include')
import numpy
# find arrayobject.h on every system an alternative would be to put
# arrayobject.h into pycogent/include, but why .. 
numpy_include_path = numpy.get_include()
distutils_extras = {"include_dirs": [include_path, numpy_include_path]}
try:
    if 'DONT_USE_PYREX' in os.environ:
        raise ImportError
    from Cython.Compiler.Version import version
    version = tuple([int(v) \
        for v in re.split("[^\d]", version) if v.isdigit()])
    if version < (0, 11, 2):
        print "Your Cython version is too old"
        raise ImportError
except ImportError:
    print "No Cython, will compile from .c files"
    for cmd in ['cython', 'pyrexc', 'predist']:
        if cmd in sys.argv:
            print "'%s' not available without Cython" % cmd
            sys.exit(1)
    from distutils.extension import Extension
    pyrex_suffix = ".c"
else:
    from Cython.Distutils import build_ext
    from Cython.Distutils.extension import Extension
    pyrex_suffix = ".pyx"
    class build_wrappers(build_ext):
        # for predist, make .c files
        def run(self):
            self.compiler = NullCompiler()
            # skip build_ext.run() and thus ccompiler setup
            build_ext.build_extensions(self)

    class build_wrappers_and_html(build_wrappers):
        def run(self):
            build_wrappers.run(self)
            build_html()

    distutils_extras["cmdclass"] = {
        'build_ext': build_ext,
        'pyrexc': build_wrappers,
        'cython': build_wrappers,
        'predist': build_wrappers_and_html}

# predist python setup.py predist --inplace --force, this is in _darcs/prefs/prefs for instructing darcs predist to execute the subsequent, predist is a darcs word

# Save some repetitive typing.  We have all compiled
# modules in place with their python siblings.
def CogentExtension(module_name, extra_compile_args=[], **kw):
    path = module_name.replace('.', '/')
    kw['extra_compile_args'] = pyrex_compile_options + extra_compile_args
    if pyrex_suffix == '.pyx':
        kw['pyrex_include_dirs'] = [include_path]
    return Extension(module_name, [path + pyrex_suffix], **kw)

short_description = "COmparative GENomics Toolkit"

# This ends up displayed by the installer
long_description = """Cogent
A toolkit for statistical analysis of biological sequences.
Version %s.
""" % __version__

setup(
    name="cogent",
    version=__version__,
    url="http://sourceforge.net/projects/pycogent",
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
    ext_modules=[
        CogentExtension("cogent.align._compare"),
        CogentExtension("cogent.align._pairwise_seqs"),
        CogentExtension("cogent.align._pairwise_pogs"),
        CogentExtension("cogent.maths._matrix_exponentiation"),
        CogentExtension("cogent.evolve._likelihood_tree"),
        CogentExtension("cogent.struct._asa"),
        CogentExtension("cogent.struct._contact"),
        CogentExtension("cogent.maths.spatial.ckd3"),

    ],
    **distutils_extras
)
