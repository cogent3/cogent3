#!/usr/bin/env python 
"""Application controller and utility functions for MSMS (molecular surface
calculation)."""
import os
import tempfile
from cogent.app.util import CommandLineApplication, ResultPath
from cogent.app.parameters import FlagParameter, ValuedParameter
from cogent.format.xyzrn import XYZRNWriter
from cogent.parse.msms import parse_VertFile


__author__ = "Marcin Cieslik"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__contributors__ = ["Marcin Cieslik"]
__license__ = "GPL"
__version__ = "1.5.3"
__maintainer__ = "Marcin Cieslik"
__email__ = "mpc4p@virginia.edu"
__status__ = "Development"


class Msms(CommandLineApplication):
    """Application controller for MSMS. The default input is a ``Entity`` 
    instance.
    
    Supported parameters:
      
      - probe_radius float : probe sphere radius, [1.5]
      - density float      : surface points density, [1.0]
      - hdensity float     : surface points high density, [3.0]
      - surface <tses,ases>: triangulated or Analytical SES, [tses]
      - no_area            : turns off the analytical surface area computation
      - noh                : ignore atoms with radius 1.2
      - no_rest_on_pbr     : no restart if pb. during triangulation
      - no_rest            : no restart if pb. are encountered
      - if filename        : sphere input file
      - of filename        : output for triangulated surface
      - af filename        : area file
      - no_header         : do not add comment line to the output
      - free_vertices      : turns on computation for isolated RS vertices
      - all_components     : compute all the surfaces components
    """

    _parameters = {
    #  -probe_radius float : probe sphere radius, [1.5]
    '-probe_radius':ValuedParameter(Prefix='-', Name='probe_radius', Delimiter=' '),
    #  -density float      : surface points density, [1.0]
    '-density':ValuedParameter(Prefix='-', Name='density', Delimiter=' '),
    #  -hdensity float     : surface points high density, [3.0]
    '-hdensity':ValuedParameter(Prefix='-', Name='hdensity', Delimiter=' '),
    #  -surface <tses,ases>: triangulated or Analytical SES, [tses]
    '-surface':ValuedParameter(Prefix='-', Name='surface', Delimiter=' '),
    #  -no_area            : turns off the analytical surface area computation
    '-no_area':FlagParameter(Prefix='-', Name='no_area', Value=False),
    #  -noh                : ignore atoms with radius 1.2      
    '-noh':FlagParameter(Prefix='-', Name='noh', Value=False),
    #  -no_rest_on_pbr     : no restart if pb. during triangulation
    '-no_rest_on_pbr':FlagParameter(Prefix='-', Name='no_rest_on_pbr', Value=False),
    #  -no_rest            : no restart if pb. are encountered
    '-no_rest':FlagParameter(Prefix='-', Name='no_rest', Value=False),
    #  -if filename        : sphere input file
    '-if':ValuedParameter(Prefix='-', Name='if', Delimiter=' ', IsPath=True),
    #  -of filename        : output for triangulated surface
    '-of':ValuedParameter(Prefix='-', Name='of', Delimiter=' ', IsPath=True),
    #  -af filename        : area file
    '-af':ValuedParameter(Prefix='-', Name='af', Delimiter=' ', IsPath=True),
    #  -no_header         : do not add comment line to the output
    '-no_header':FlagParameter(Prefix='-', Name='no_header', Value=False),
    #  -free_vertices      : turns on computation for isolated RS vertices
    '-free_vertices':FlagParameter(Prefix='-', Name='free_vertices', Value=False),
    #  -all_components     : compute all the surfaces components
    '-all_components':FlagParameter(Prefix='-', Name='all_components', Value=False),
    #######################
    #  -one_cavity #atoms at1 [at2][at3] : Compute the surface for an internal  
    # cavity for which at least one atom is specified
    #######################               
    #  -socketName servicename : socket connection from a client
    #  -socketPort portNumber : socket connection from a client
    #  -xdr                : use xdr encoding over socket
    #  -sinetd             : inetd server connectio
    }

    _command = "msms"
    _input_handler = '_input_from_entity'

    def _input_from_entity(self, data):
        """This allows to feed entities to msms."""
        if data:
            # create temporary files and names.
            fd, self._input_filename = tempfile.mkstemp()
            os.close(fd)
            # write XYZR data
            fh = open(self._input_filename, 'wb')
            XYZRNWriter(fh, data)
            fh.close()
            #
            self.Parameters['-if'].on(self._input_filename)
            self.Parameters['-of'].on(self._input_filename) # msms appends .vert .face
            self.Parameters['-af'].on(self._input_filename) # msms appends .area
        if (not self.Parameters['-if'].isOn()) or \
           (not self.Parameters['-of'].isOn()) or \
           (not self.Parameters['-af'].isOn()):
            raise ValueError('All of -if, -of and -af have to be specified.')
        return ""

    def _get_result_paths(self, data):
        result = {}
        vert_file = self.Parameters['-of'].Value + '.vert'
        result['VertFile'] = ResultPath(Path=vert_file, IsWritten=True)
        face_file = self.Parameters['-of'].Value + '.face'
        result['FaceFile'] = ResultPath(Path=face_file, IsWritten=True)
        if not self.Parameters['-no_area'].Value:
            area_file = self.Parameters['-af'].Value + '.area'
            result['AreaFile'] = ResultPath(Path=area_file, IsWritten=True)
        return result
    
    
def surface_xtra(entity, **kwargs):
    """Runs the MSMS application to create the molecular surface, which is an
    a numpy array of 3D-coordinates.
    
    Arguments:
    
        * entity - an ``Entity`` instance.
    
    Additional keyworded arguments are for the ``Msms`` application controller.
    Returns the numpy array and put it into entity.xtra['MSMS_SURFACE'].  
    """
    msms = Msms(**kwargs)
    res = msms(entity)
    surface = parse_VertFile(res['VertFile'])
    entity.xtra['MSMS_SURFACE'] = surface
    res.cleanUp()   # remove all temporary files
    return surface
