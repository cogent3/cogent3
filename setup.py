#!/usr/bin/env python
from distutils.core import setup
from distutils.extension import Extension
import sys, os, re

__author__ = "Peter Maxwell"
__copyright__ = "Copyright 2007-2008, The Cogent Project"
__contributors__ = ["Peter Maxwell", "Gavin Huttley", "Matthew Wakefield", 
                    "Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.1"
__maintainer__ = "Peter Maxwell"
__email__ = "pm67nz@gmail.com"
__status__ = "Production"

doc_imports_failed = False
try:
    import docutils.core
    from docpytils.docpytils import *
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
        print "Failed to build html due to ImportErrors for either docutils or docpytils"
        return
    stylesheet_path = os.path.join(os.getcwd(), "doc", "html_style.css")
    new_html = docutils.core.publish_file(source_path="doc/index.rest",
                    writer_name='html',
                    destination_path='doc/index.html',
                    settings_overrides={"embed-stylesheet": True,
                                        "stylesheet_path": stylesheet_path})
    print "Built index.html"

# Compiling Pyrex modules to .c and .so
try:
    if 'DONT_USE_PYREX' in os.environ:
        raise ImportError
    else:
        from Pyrex.Compiler.Version import version
        version = tuple([int(v) \
            for v in re.split("[^\d]",version) if v.isdigit()])
        if version < (0, 9, 4):
            print "Your Pyrex version is too old"
            raise ImportError
except ImportError:
    # build from intermediate .c files
    # if we don't have pyrex
    print "Didn't find Pyrex - will compile from .c files"
    distutils_extras = {"include_dirs": [os.path.join(os.getcwd(), 'include')]}
    pyrex_suffix = ".c"
else:
    from Pyrex.Distutils import build_ext
                        
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

    distutils_extras = {
        "cmdclass": {
            'build_ext': build_ext,
            'pyrexc': build_wrappers,
            'predist': build_wrappers_and_html},
        "include_dirs": [os.path.join(os.getcwd(), 'include')],
        }
    
    pyrex_suffix = ".pyx"


# predist python setup.py predist --inplace --force, this is in _darcs/prefs/prefs for instructing darcs predist to execute the subsequent, predist is a darcs word

# Save some repetitive typing.  We have all compiled
# modules in place with their python siblings.
def CogentExtension(module_name, extra_compile_args=[], include_dirs=[], **kw):
    path = module_name.replace('.', '/')
    kw['extra_compile_args'] = pyrex_compile_options + extra_compile_args
    return Extension(module_name, [path + pyrex_suffix], **kw)

# get version
execfile('cogent/__init__.py')

short_description = "COmparative GENomics Toolkit"

# This ends up displayed by the installer
long_description = """Cogent
A toolkit for statistical analysis of biological sequences.
Version %s.
""" % version

setup(
    name = "cogent",
    version = version,
    url = "http://sourceforge.net/projects/pycogent",
    author = "Gavin Huttley, Rob Knight",
    author_email = "gavin.huttley@anu.edu.au, rob@spot.colorado.edu",
    description = short_description,
    long_description = long_description,
    platforms = ["any"],
    license = ["GPL"],
    keywords = ["biology", "genomics", "statistics", "phylogeny", "evolution",
                "bioinformatics"],
    classifiers = [
            "Development Status :: 5 - Production/Stable",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: GNU General Public License (GPL)",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            "Topic :: Software Development :: Libraries :: Python Modules",
            "Operating System :: OS Independent",
            ],
    packages = ['cogent', 'cogent.align', 'cogent.align.weights', 'cogent.app',
                'cogent.cluster', 'cogent.core', 'cogent.data', 'cogent.db',
                'cogent.draw', 'cogent.draw.matplotlib', 'cogent.evolve',
                'cogent.format', 'cogent.maths', 'cogent.maths.matrix',
                'cogent.maths.stats', 'cogent.maths.stats.cai', 'cogent.motif',
                'cogent.parse', 'cogent.phylo', 'cogent.recalculation',
                'cogent.seqsim', 'cogent.struct', 'cogent.util'],
    ext_modules=[
        CogentExtension("cogent.align._compare"),
        CogentExtension("cogent.align._pairwise_seqs"),
        CogentExtension("cogent.align._pairwise_pogs"),
        CogentExtension("cogent.maths._matrix_exponentiation"),
        CogentExtension("cogent.evolve._likelihood_tree"),
    ],
    **distutils_extras
)
