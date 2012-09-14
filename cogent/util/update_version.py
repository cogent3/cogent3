#!/usr/bin/env python

"""Support for updating version strings in the PyCogent source tree

All .py, .pyx, and .c files descending from cogent/ will be updated
All .py files descending from tests/ will be updated
All .pyx and .h files descending from include/ will be updated
docs/conf.py will be updated
setup.py will be updated
###cogent-requirements.txt will be updated### not being updated now
"""

from optparse import make_option, OptionParser
from sys import argv
from os import path
import os

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

options = [make_option('--pycogent_dir',dest='pycogent_dir',type='string',
                       default=''),
           make_option('--new_version',dest='version',type='string',
                       default=''),
           make_option('--is_release',dest='is_release',\
                       action='store_true', default=False),
           make_option('--verbose',dest='verbose',action='store_true',
                       default=False),
           make_option('--mock_run',dest='mockrun',action='store_true',
                       default=False),
           make_option('--new_version_short',dest='version_short',type='string',
                       default=None)]

class VersionUpdater(object):
    """Handles version update of files contained within the PyCogent tree"""
    def __init__(self, PyCogentDirectory=None, Version=None, \
            IsRelease=False, Verbose=False, MockRun=False, VersionShort=None):
        self.PyCogentDirectory = PyCogentDirectory
        self.Version = Version
        self.VersionShort = VersionShort
        self.VersionTuple = tuple(self.Version.split('.'))
        self.IsRelease = IsRelease
        self.Verbose = Verbose
        self.MockRun = MockRun
        
        self.CodesDirectory = path.join(self.PyCogentDirectory, 'cogent')
        self.TestsDirectory = path.join(self.PyCogentDirectory, 'tests')
        self.DocDirectory = path.join(self.PyCogentDirectory, 'doc')
        self.IncludesDirectory = path.join(self.PyCogentDirectory, 'include')
        
        if not os.access(path.join(self.CodesDirectory, '__init__.py'),os.R_OK):
            raise IOError, "Could not locate cogent/__init__.py"
        if not os.access(path.join(self.TestsDirectory, '__init__.py'),os.R_OK):
            raise IOError, "Could not locate tests/__init__.py"
        if not os.access(path.join(self.DocDirectory, 'conf.py'), os.R_OK):
            raise IOError, "Could not locate doc/conf.py"
        if not os.access(path.join(self.IncludesDirectory, \
                'array_interface.h'), os.R_OK):
            raise IOError, "Cound not locate include/array_interface.h"

    def _get_base_files(self):
        """Support method, provides relative locations for files in base dir"""
        setup_file = path.join(self.PyCogentDirectory, 'setup.py')
        #reqs_file = path.join(self.PyCogentDirectory, 'cogent-requirements.txt')
        #return [(setup_file, 'Python'), (reqs_file, 'Properties')]
        return [(setup_file, 'Python')]

    def _get_test_files(self):
        """Support method, provides relative locations for test files"""
        for dirpath, dirnames, filenames in os.walk(self.TestsDirectory):
            for f in filenames:
                if f.endswith('.py'):
                    yield (path.join(dirpath, f), 'Python')

    def _get_code_files(self):
        """Support method, provides relative locations for code files
        
        Yields file name and file type
        """
        for dirpath, dirnames, filenames in os.walk(self.CodesDirectory):
            for f in filenames:
                rel_name = path.join(dirpath, f)
                if f.endswith('.py'):
                    yield (rel_name, 'Python')
                elif f.endswith('.pyx'):
                    yield (rel_name, 'PyRex')
                elif f.endswith('.c'):
                    yield (rel_name, 'C')
                else:
                    pass

    def _get_doc_files(self):
        """Support method, provides relative locations for test files
        
        Only yields conf.py currently
        """
        return [(path.join(self.DocDirectory, 'conf.py'), 'Python')]

    def _get_include_files(self):
        """Support method, provides relative locations for include files

        Yields file name and file type
        """
        for dirpath, dirnames, filenames in os.walk(self.IncludesDirectory):
            for f in filenames:
                rel_name = path.join(dirpath, f)
                if f.endswith('.pyx'):
                    yield (rel_name, 'PyRex')
                elif f.endswith('.h'):
                    yield (rel_name, 'Header')
                else:
                    pass

    def _update_python_file(self, lines, filename):
        """Updates the __version__ string of a Python file"""
        found_version_line = False
        for lineno, line in enumerate(lines):
            if line.startswith('__version__'):
                found_version_line = True
                break
        if found_version_line:
            if self.Verbose:
                print 'Version string found on line %d' % lineno
            lines[lineno] = '__version__ = "%s"\n' % self.Version
        else:
            print "No version string found in %s" % filename
        return (lines, found_version_line)

    def _update_properties_file(self, lines, filename):
        """Updates version information in specific properties files

        Expects the properties file to be in "key=value" lines
        """
        found_version_line = False
        if filename.endswith('cogent-requirements.txt'):
            for lineno, line in enumerate(lines):
                if 'packages/source/c/cogent' in line:
                    found_version_line = True
                    break
        if found_version_line:
            if self.Verbose:
                print 'Version string found on line %d' % lineno
            http_base = lines[lineno].rsplit('/',1)[0]
            lines[lineno] = '%s/PyCogent-%s.tgz\n' % (http_base, self.Version)
        else:
            print "No version string found in %s" % filename
        return (lines, found_version_line) 

    def _update_doc_conf_file(self, lines, filename):
        """Updates doc/conf.py file"""
        versionline = None
        releaseline = None

        for lineno, line in enumerate(lines):
            if line.startswith('version'):
                versionline = lineno
            if line.startswith('release'):
                releaseline = lineno
            if versionline is not None and releaseline is not None:
                break

        if versionline is None:
            print "No version string found in doc/conf.py"
        else:
            if self.Verbose:
                print 'Version string found on line %d' % versionline
            lines[versionline] = 'version = "%s"\n' % self.VersionShort

        if releaseline is None:
            print "No release string found in doc/conf.py"
        else:
            if self.Verbose:
                print 'Release string found on line %d' % releaseline
            lines[releaseline] = 'release = "%s"\n' % self.Version

        return (lines, versionline and releaseline)

    def _update_pyrex_file(self, lines, filename):
        """Updates __version__ within a pyx file"""
        found_version_line = False
        for lineno, line in enumerate(lines):
            if line.startswith('__version__'):
                found_version_line = True
                break
        if found_version_line:
            if self.Verbose:
                print 'Version string found on line %d' % lineno
            lines[lineno] = '__version__ = "%s"\n' % str(self.VersionTuple)
        else:
            print "No version string found in %s" % filename
        return (lines, found_version_line)

    def _update_header_file(self, lines, filename):
        """Updates a C header file"""
        found_version_line = False
        for lineno, line in enumerate(lines):
            if line.startswith('#define PYCOGENT_VERSION'):
                found_version_line = True
                break
        if found_version_line:
            if self.Verbose:
                print 'Version string found on line %d' % lineno
            lines[lineno] = '#define PYCOGENT_VERSION "%s"\n' \
                    % self.Version
        else:
            print "No version string found in %s" % filename
        return (lines, found_version_line)

    def _update_c_file(self, lines, filename):
        """Updates a C file"""
        # same as C header...
        return self._update_header_file(lines, filename)

    def _file_writer(self, lines, filename):
        """Handles writing out to the file system"""
        if self.MockRun:
            return

        if self.Verbose:
            print "Writing file %s" % filename

        updated_file = open(filename, 'w')
        updated_file.write(''.join(lines))
        updated_file.close()

    def updateBaseFiles(self):
        """Updates version strings in files in base PyCogent directory"""
        for filename, filetype in self._get_base_files():
            lines = open(filename).readlines()

            if self.Verbose:
                print 'Reading %s' % filename

            if filetype is 'Python':
                lines, write_out = self._update_python_file(lines, filename) 
            elif filetype is 'Properties':
                lines, write_out = self._update_properties_file(lines,filename)
            else:
                raise TypeError, "Unknown base file type %s" % filetype

            if write_out:
                self._file_writer(lines, filename) 

    def updateDocFiles(self):
        """Updates version strings in documentation files
        
        So far we only update conf.py
        """
        for filename, filetype in self._get_doc_files():
            lines = open(filename).readlines()

            if self.Verbose:
                print 'Reading %s' % filename
           
            if filename.endswith('conf.py'):
                lines, write_out = self._update_doc_conf_file(lines, filename)
            else:
                raise TypeError, "Unknown doc file type: %s" % filetype

            if write_out:
                self._file_writer(lines, filename) 
        
    def updateIncludeFiles(self):
        """Updates version strings in include files"""
        for filename, filetype in self._get_include_files():
            lines = open(filename).readlines()
            found_version_line = False

            if self.Verbose:
                print 'Reading %s' % filename
            
            if filetype is 'PyRex':
                lines, write_out = self._update_pyrex_file(lines, filename)
            elif filetype is 'Header':
                lines, write_out = self._update_header_file(lines, filename)
            else:
                raise TypeError, "Unknown include file type %s" % filetype

            if write_out:
                self._file_writer(lines, filename) 

    def updateTestFiles(self):
        """Updates version strings in test files"""
        for filename, filetype in self._get_test_files():
            lines = open(filename).readlines()
            found_version_line = False

            if self.Verbose:
                print 'Reading %s' % filename

            if filetype is 'Python':
                lines, write_out = self._update_python_file(lines, filename)
            else:
                raise TypeError, "Unknown test file type %s" % filetype

            if write_out:
                self._file_writer(lines, filename) 

    def updateCodeFiles(self):
        """Updates version strings in code files"""
        # if this annoying slow, could probably drop to bash or soemthing
        # for a search/replace
        for filename, filetype in self._get_code_files():
            lines = open(filename).readlines()
            found_version_line = False

            if self.Verbose:
                print 'Reading %s' % filename

            if filetype is 'Python':
                lines, write_out = self._update_python_file(lines, filename)
            elif filetype is 'PyRex':
                lines, write_out = self._update_pyrex_file(lines, filename)
            elif filetype is 'C':
                lines, write_out = self._update_c_file(lines, filename)
            else:
                raise TypeError, "Unknown code file type %s" % filetype

            if write_out:
                self._file_writer(lines, filename) 

def main(arg_list=argv):
    parser = OptionParser(option_list=options)
    opts, args = parser.parse_args(args=arg_list)

    updater = VersionUpdater(PyCogentDirectory=opts.pycogent_dir,
                             Version=opts.version,
                             VersionShort=opts.version_short,
                             IsRelease=opts.is_release,
                             Verbose=opts.verbose,
                             MockRun=opts.mockrun)

    updater.updateCodeFiles()
    updater.updateTestFiles()
    updater.updateDocFiles()
    updater.updateIncludeFiles()
    updater.updateBaseFiles()

if __name__ == '__main__':
    main(argv)
