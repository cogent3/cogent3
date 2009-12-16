#!/usr/bin/env python

from cogent.util.unit_test import TestCase, main
from cogent.app.rnaview import RnaView
from tempfile import mktemp, tempdir
from os import remove, system, getcwd

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Greg Caporaso", "Sandra Smit", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.4"
__maintainer__ = "Greg Caporaso"
__email__ = "Greg Caporaso"
__status__ = "Production"

class test_rnaview(TestCase):
    """ Tests of the rnaview application controller """

    def setUp(self):
        """ SetUp some objects for use by tests """
        # Define the directory where mktemp creates files, this simplifies 
        # things b/c if ~/tmp exists files are created there, otherwise they 
        # are created in /tmp. For the case of tests, it's easier just to 
        # set it to a constant
        self.fake_pdb_file = fake_pdb_file
        self.fake_nmr_file = fake_nmr_file
        self.r1 = RnaView(WorkingDir='/tmp')
        self.r_input_as_lines = RnaView(WorkingDir='/tmp',\
            InputHandler='_input_as_lines')
        self.r2 = RnaView(params={'-v':None},WorkingDir='/tmp')
        self.r3 = RnaView(params={'-c':'A'},WorkingDir='/tmp')
        self.r4 = RnaView(params={'-x':None,'-p':None},WorkingDir='/tmp')
        self.r5 = RnaView()

    def test_base_command(self):
        """RnaView: BaseCommand is built correctly """
        self.assertEqual(self.r1.BaseCommand, 'cd "/tmp/"; rnaview')
        self.assertEqual(self.r2.BaseCommand, 'cd "/tmp/"; rnaview -v')
        self.assertEqual(self.r3.BaseCommand, 'cd "/tmp/"; rnaview -c A')
        assert (self.r4.BaseCommand == 'cd "/tmp/"; rnaview -x -p') or \
            (self.r4.BaseCommand == 'cd "/tmp/"; rnaview -p -x')
        self.assertEqual(self.r5.BaseCommand,\
            'cd "' + getcwd() + '/'+ '"; rnaview')

    def test_file_pointers_no_extras_input_as_file(self):
        """RnaView: pointers created only for minimal files w/ _input_as_lines
        """
        written_files = {}.fromkeys(['bp_stats','base_pairs',\
            'StdOut','StdErr','ExitStatus'])

        res =self.r_input_as_lines(data=self.fake_pdb_file.split('\n'))
        for f in res:
            if f in written_files:
                assert res[f] is not None
            else:
                assert res[f] is None
        res.cleanUp()
 
    def test_file_pointers_no_extras(self):
        """RnaView: pointers created only for minimal files
        """
        filename = mktemp() + '.pdb'
        f = open(filename,"w")
        f.writelines(self.fake_pdb_file)
        f.close()
        
        written_files = {}.fromkeys(['bp_stats','base_pairs',\
            'StdOut','StdErr','ExitStatus'])
        res = self.r1(data=filename)
        for f in res:
            if f in written_files:
                assert res[f] is not None
            else:
                assert res[f] is None
        res.cleanUp()
        remove(filename)
   
    def test_file_pointers_w_vrml(self):
        """RnaView: pointers created for minimal files and files turned on
        """
        # need to make a fake pdb file with base pairs so that wrl file will be 
        # created
        filename = mktemp() + '.pdb'
        f = open(filename,"w")
        f.writelines(self.fake_pdb_file)
        f.close()

        written_files = {}.fromkeys(['bp_stats','base_pairs',\
            'StdOut','StdErr','ExitStatus','vrml'])
        res = self.r2(data=filename)
        for f in res:
            if f in written_files:
                assert res[f] is not None
            else:
                assert res[f] is None
        
        res.cleanUp()
        remove(filename)

    def test_base_pairs_out(self):
        """RnaView: output sanity check """
        filename = mktemp()
        f = open(filename,"w")
        f.writelines(self.fake_pdb_file)
        f.close()
        
        res = self.r2(data=filename)
        bp_file = list(res['base_pairs'])
        self.assertEqual(bp_file[0],'PDB data file name: ' + filename + '\n')
        
        res.cleanUp()
        remove(filename)
        
    def test_get_pdb_filename(self):
        """RnaView: _get_pdb_filename functions as expected """
        #Test X-RAY CRYSTALLOGRAPHY file
        f = open('/tmp/1EHZ.pdb',"w")
        f.writelines(self.fake_pdb_file)
        f.close()
        self.assertEqual(self.r1._get_pdb_filename('/tmp/1EHZ.pdb'),'1EHZ.pdb')
        remove('/tmp/1EHZ.pdb')
        
        #Test NMR file
        f = open('/tmp/pdb17.ent',"w")
        f.writelines(self.fake_nmr_file)
        f.close()
        self.assertEqual(self.r1._get_pdb_filename('/tmp/pdb17.ent'),\
            'pdb17.ent_nmr.pdb')
        remove('/tmp/pdb17.ent')
        
        #OLD TESTS BEFORE THE ACTUAL FILE HAD TO BE THERE
        #self.assertEqual(self.r1._get_pdb_filename('/tmp/1EHZ.pdb'),'1EHZ.pdb')
        #self.assertEqual(self.r1._get_pdb_filename('/tmp/duh/1EHZ.pdb'),\
        #    '1EHZ.pdb')
        #self.assertEqual(self.r1._get_pdb_filename('1EHZ1.pdb'),'1EHZ1.pdb')
        #self.assertEqual(self.r1._get_pdb_filename('/1'),'1')
        #self.assertEqual(self.r1._get_pdb_filename('1EHZZZ.pdb'),'1EHZZZ.pdb')
        #self.assertEqual(self.r1._get_pdb_filename('/tmp/tmpW3urtc'),\
        #    'tmpW3urtc')

    def test_get_out_path(self):
        """RnaView: _get_out_path functions as expected """
        self.assertEqual(self.r1._get_out_path('1EHZ.pdb'),'')
        self.assertEqual(self.r1._get_out_path('/tmp/1EHZ.pdb'),'/tmp/')
        self.assertEqual(self.r1._get_out_path('/tmp/duh/1EHZ.pdb'),'/tmp/duh/')
        self.assertEqual(self.r1._get_out_path('/1'),'/')


    def test_accept_exit_status(self):
        """RnaView: _accept_exit_status functions as expected """
        self.assertEqual(self.r1._accept_exit_status(0),True)
        self.assertEqual(self.r1._accept_exit_status(1),False)
        self.assertEqual(self.r1._accept_exit_status('0'),False)
        self.assertEqual(self.r1._accept_exit_status(None),False)
        self.assertEqual(self.r1._accept_exit_status(''),False)

   

fake_pdb_file = """ATOM      1  O3P   G A   1      50.193  51.190  50.534  1.00 99.85           O  
ATOM      2  P     G A   1      50.626  49.730  50.573  1.00100.19           P  
ATOM      3  O1P   G A   1      49.854  48.893  49.562  1.00100.19           O  
ATOM      4  O2P   G A   1      52.137  49.542  50.511  1.00 99.21           O  
ATOM      5  O5*   G A   1      50.161  49.136  52.023  1.00 99.82           O  
ATOM      6  C5*   G A   1      50.216  49.948  53.210  1.00 98.63           C  
ATOM      7  C4*   G A   1      50.968  49.231  54.309  1.00 97.84           C  
ATOM      8  O4*   G A   1      50.450  47.888  54.472  1.00 97.10           O  
ATOM      9  C3*   G A   1      52.454  49.030  54.074  1.00 98.07           C  
ATOM     10  O3*   G A   1      53.203  50.177  54.425  1.00 99.39           O  
ATOM     11  C2*   G A   1      52.781  47.831  54.957  1.00 96.96           C  
ATOM     12  O2*   G A   1      53.018  48.156  56.313  1.00 96.77           O  
ATOM     13  C1*   G A   1      51.502  47.007  54.836  1.00 95.70           C  
ATOM     14  N9    G A   1      51.628  45.992  53.798  1.00 93.67           N  
ATOM     15  C8    G A   1      51.064  46.007  52.547  1.00 92.60           C  
ATOM     16  N7    G A   1      51.379  44.966  51.831  1.00 91.19           N  
ATOM     17  C5    G A   1      52.197  44.218  52.658  1.00 91.47           C  
ATOM     18  C6    G A   1      52.848  42.992  52.425  1.00 90.68           C  
ATOM     19  O6    G A   1      52.826  42.291  51.404  1.00 90.38           O  
ATOM     20  N1    G A   1      53.588  42.588  53.534  1.00 90.71           N  
ATOM     21  C2    G A   1      53.685  43.282  54.716  1.00 91.21           C  
ATOM     22  N2    G A   1      54.452  42.733  55.671  1.00 91.23           N  
ATOM     23  N3    G A   1      53.077  44.429  54.946  1.00 91.92           N  
ATOM     24  C4    G A   1      52.356  44.836  53.879  1.00 92.62           C  
ATOM     25  P     C A   2      54.635  50.420  53.741  1.00100.19           P  
ATOM     26  O1P   C A   2      55.145  51.726  54.238  1.00100.19           O  
ATOM     27  O2P   C A   2      54.465  50.204  52.269  1.00100.19           O  
ATOM     28  O5*   C A   2      55.563  49.261  54.342  1.00 98.27           O  
ATOM     29  C5*   C A   2      55.925  49.246  55.742  1.00 95.40           C  
ATOM     30  C4*   C A   2      56.836  48.075  56.049  1.00 93.33           C  
ATOM     31  O4*   C A   2      56.122  46.828  55.830  1.00 92.18           O  
ATOM     32  C3*   C A   2      58.090  47.947  55.197  1.00 92.75           C  
ATOM     33  O3*   C A   2      59.174  48.753  55.651  1.00 92.89           O  
ATOM     34  C2*   C A   2      58.416  46.463  55.298  1.00 91.81           C  
ATOM     35  O2*   C A   2      59.140  46.136  56.466  1.00 91.36           O  
ATOM     36  C1*   C A   2      57.022  45.836  55.356  1.00 90.59           C  
ATOM     37  N1    C A   2      56.570  45.364  54.029  1.00 88.84           N  
ATOM     38  C2    C A   2      57.094  44.157  53.520  1.00 88.64           C  
ATOM     39  O2    C A   2      57.921  43.516  54.198  1.00 88.97           O  
ATOM     40  N3    C A   2      56.686  43.721  52.301  1.00 87.36           N  
ATOM     41  C4    C A   2      55.802  44.437  51.597  1.00 87.11           C  
ATOM     42  N4    C A   2      55.430  43.972  50.397  1.00 86.30           N  
ATOM     43  C5    C A   2      55.259  45.660  52.089  1.00 86.87           C  
ATOM     44  C6    C A   2      55.663  46.080  53.296  1.00 88.01           C  
ATOM     45  P     G A   3      60.184  49.419  54.574  1.00 92.31           P  
ATOM     46  O1P   G A   3      61.015  50.422  55.295  1.00 92.97           O  
ATOM     47  O2P   G A   3      59.371  49.857  53.404  1.00 91.56           O  
ATOM     48  O5*   G A   3      61.137  48.219  54.105  1.00 88.57           O  
ATOM     49  C5*   G A   3      62.175  47.724  54.969  1.00 83.44           C  
ATOM     50  C4*   G A   3      62.769  46.443  54.422  1.00 79.87           C  
ATOM     51  O4*   G A   3      61.734  45.427  54.299  1.00 78.36           O  
ATOM     52  C3*   G A   3      63.405  46.499  53.040  1.00 78.97           C  
ATOM     53  O3*   G A   3      64.741  47.029  53.060  1.00 79.76           O  
ATOM     54  C2*   G A   3      63.359  45.032  52.608  1.00 77.19           C  
ATOM     55  O2*   G A   3      64.411  44.256  53.155  1.00 77.80           O  
ATOM     56  C1*   G A   3      62.018  44.572  53.194  1.00 73.98           C  
ATOM     57  N9    G A   3      60.934  44.675  52.202  1.00 68.20           N  
ATOM     58  C8    G A   3      60.024  45.702  52.050  1.00 65.03           C  
ATOM     59  N7    G A   3      59.252  45.556  51.003  1.00 62.99           N  
ATOM     60  C5    G A   3      59.655  44.348  50.447  1.00 59.95           C  
ATOM     61  C6    G A   3      59.189  43.675  49.292  1.00 55.65           C  
ATOM     62  O6    G A   3      58.287  44.013  48.522  1.00 53.32           O  
ATOM     63  N1    G A   3      59.893  42.491  49.072  1.00 54.00           N  
ATOM     64  C2    G A   3      60.906  42.006  49.876  1.00 55.46           C  
ATOM     65  N2    G A   3      61.512  40.873  49.479  1.00 48.16           N  
ATOM     66  N3    G A   3      61.312  42.605  50.983  1.00 56.69           N  
ATOM     67  C4    G A   3      60.666  43.774  51.193  1.00 61.76           C  
ATOM     68  P     G A   4      65.295  47.868  51.793  1.00 79.34           P  
ATOM     69  O1P   G A   4      66.538  48.562  52.246  1.00 80.87           O  
ATOM     70  O2P   G A   4      64.193  48.679  51.209  1.00 79.00           O  
ATOM     71  O5*   G A   4      65.720  46.752  50.724  1.00 75.17           O  
ATOM     72  C5*   G A   4      66.789  45.843  51.019  1.00 68.95           C  
ATOM     73  C4*   G A   4      66.749  44.634  50.114  1.00 65.13           C  
ATOM     74  O4*   G A   4      65.484  43.939  50.258  1.00 61.83           O  
ATOM     75  C3*   G A   4      66.881  44.840  48.611  1.00 62.79           C  
ATOM     76  O3*   G A   4      68.230  44.977  48.176  1.00 61.75           O  
ATOM     77  C2*   G A   4      66.318  43.538  48.064  1.00 60.58           C  
ATOM     78  O2*   G A   4      67.283  42.514  48.122  1.00 59.59           O  
ATOM     79  C1*   G A   4      65.192  43.241  49.051  1.00 58.29           C  
ATOM     80  N9    G A   4      63.923  43.716  48.500  1.00 53.36           N  
ATOM     81  C8    G A   4      63.204  44.843  48.842  1.00 49.19           C  
ATOM     82  N7    G A   4      62.140  45.013  48.107  1.00 46.88           N  
ATOM     83  C5    G A   4      62.144  43.926  47.243  1.00 45.95           C  
ATOM     84  C6    G A   4      61.246  43.573  46.206  1.00 43.95           C  
ATOM     85  O6    G A   4      60.182  44.136  45.874  1.00 42.46           O  
ATOM     86  N1    G A   4      61.672  42.428  45.540  1.00 40.34           N  
ATOM     87  C2    G A   4      62.788  41.686  45.867  1.00 42.60           C  
ATOM     88  N2    G A   4      63.034  40.588  45.135  1.00 38.91           N  
ATOM     89  N3    G A   4      63.612  41.994  46.850  1.00 44.40           N  
ATOM     90  C4    G A   4      63.239  43.117  47.480  1.00 48.65           C  
ATOM     91  P     A A   5      68.530  45.722  46.789  1.00 59.03           P  
ATOM     92  O1P   A A   5      69.991  45.842  46.548  1.00 60.84           O  
ATOM     93  O2P   A A   5      67.685  46.959  46.834  1.00 60.64           O  
ATOM     94  O5*   A A   5      67.957  44.735  45.675  1.00 57.14           O  
ATOM     95  C5*   A A   5      68.648  43.529  45.323  1.00 53.41           C  
ATOM     96  C4*   A A   5      67.927  42.844  44.191  1.00 50.63           C  
ATOM     97  O4*   A A   5      66.589  42.480  44.646  1.00 48.70           O  
ATOM     98  C3*   A A   5      67.665  43.715  42.964  1.00 50.77           C  
ATOM     99  O3*   A A   5      68.747  43.769  42.051  1.00 52.86           O  
ATOM    100  C2*   A A   5      66.455  43.024  42.355  1.00 48.94           C  
ATOM    101  O2*   A A   5      66.864  41.798  41.731  1.00 48.54           O  
ATOM    102  C1*   A A   5      65.646  42.719  43.615  1.00 44.50           C  
ATOM    103  N9    A A   5      64.779  43.843  44.021  1.00 42.01           N  
ATOM    104  C8    A A   5      64.938  44.803  45.016  1.00 39.75           C  
ATOM    105  N7    A A   5      63.925  45.649  45.113  1.00 41.58           N  
ATOM    106  C5    A A   5      63.049  45.220  44.115  1.00 38.26           C  
ATOM    107  C6    A A   5      61.796  45.688  43.683  1.00 35.83           C  
ATOM    108  N6    A A   5      61.110  46.688  44.232  1.00 32.66           N  
ATOM    109  N1    A A   5      61.233  45.057  42.644  1.00 35.14           N  
ATOM    110  C2    A A   5      61.870  44.017  42.074  1.00 38.97           C  
ATOM    111  N3    A A   5      63.024  43.467  42.399  1.00 36.02           N  
ATOM    112  C4    A A   5      63.571  44.119  43.437  1.00 39.04           C  
ATOM    113  P     U A   6      69.150  45.179  41.392  1.00 55.09           P  
ATOM    114  O1P   U A   6      70.511  44.926  40.836  1.00 56.37           O  
ATOM    115  O2P   U A   6      68.953  46.283  42.381  1.00 51.00           O  
ATOM    116  O5*   U A   6      68.119  45.358  40.184  1.00 50.38           O  
ATOM    117  C5*   U A   6      67.912  44.271  39.258  1.00 48.10           C  
ATOM    118  C4*   U A   6      66.579  44.400  38.565  1.00 47.17           C  
ATOM    119  O4*   U A   6      65.486  44.324  39.513  1.00 46.61           O  
ATOM    120  C3*   U A   6      66.344  45.708  37.850  1.00 45.27           C  
ATOM    121  O3*   U A   6      66.964  45.696  36.590  1.00 45.77           O  
ATOM    122  C2*   U A   6      64.833  45.733  37.727  1.00 45.88           C  
ATOM    123  O2*   U A   6      64.431  44.864  36.684  1.00 44.33           O  
ATOM    124  C1*   U A   6      64.413  45.113  39.057  1.00 41.32           C  
ATOM    125  N1    U A   6      64.065  46.111  40.079  1.00 39.88           N  
ATOM    126  C2    U A   6      62.798  46.658  39.977  1.00 36.06           C  
ATOM    127  O2    U A   6      62.021  46.333  39.099  1.00 38.25           O  
ATOM    128  N3    U A   6      62.487  47.582  40.924  1.00 34.15           N  
ATOM    129  C4    U A   6      63.272  48.002  41.975  1.00 37.36           C  
ATOM    130  O4    U A   6      62.822  48.829  42.752  1.00 39.30           O  
ATOM    131  C5    U A   6      64.583  47.395  42.032  1.00 39.23           C  
ATOM    132  C6    U A   6      64.926  46.497  41.084  1.00 35.72           C  
ATOM    133  P     U A   7      67.463  47.074  35.969  1.00 44.37           P  
ATOM    134  O1P   U A   7      68.318  46.756  34.822  1.00 48.09           O  
ATOM    135  O2P   U A   7      67.945  47.948  37.077  1.00 45.68           O  
ATOM    136  O5*   U A   7      66.104  47.724  35.455  1.00 40.88           O  
ATOM    137  C5*   U A   7      65.285  47.024  34.459  1.00 37.89           C  
ATOM    138  C4*   U A   7      64.055  47.852  34.101  1.00 35.74           C  
ATOM    139  O4*   U A   7      63.297  48.107  35.326  1.00 38.13           O  
ATOM    140  C3*   U A   7      64.317  49.197  33.459  1.00 36.87           C  
ATOM    141  O3*   U A   7      63.402  49.394  32.378  1.00 37.45           O  
ATOM    142  C2*   U A   7      64.097  50.171  34.624  1.00 36.55           C  
ATOM    143  O2*   U A   7      63.595  51.417  34.246  1.00 35.54           O  
ATOM    144  C1*   U A   7      63.015  49.475  35.442  1.00 37.23           C  
ATOM    145  N1    U A   7      63.056  49.858  36.864  1.00 36.91           N  
ATOM    146  C2    U A   7      62.011  50.628  37.343  1.00 34.52           C  
ATOM    147  O2    U A   7      61.087  50.966  36.653  1.00 34.65           O  
ATOM    148  N3    U A   7      62.112  50.993  38.659  1.00 37.03           N  
ATOM    149  C4    U A   7      63.131  50.684  39.541  1.00 40.15           C  
ATOM    150  O4    U A   7      63.105  51.143  40.699  1.00 36.62           O  
ATOM    151  C5    U A   7      64.179  49.865  38.971  1.00 36.52           C  
ATOM    152  C6    U A   7      64.106  49.490  37.691  1.00 36.25           C  
ATOM    153  P     U A   8      63.884  49.282  30.858  1.00 36.77           P  
ATOM    154  O1P   U A   8      62.852  49.899  29.952  1.00 38.95           O  
ATOM    155  O2P   U A   8      64.442  47.955  30.547  1.00 38.70           O  
ATOM    156  O5*   U A   8      65.171  50.254  30.733  1.00 35.95           O  
ATOM    157  C5*   U A   8      64.994  51.676  30.500  1.00 33.53           C  
ATOM    158  C4*   U A   8      66.105  52.236  29.628  1.00 34.33           C  
ATOM    159  O4*   U A   8      67.428  52.119  30.261  1.00 31.81           O  
ATOM    160  C3*   U A   8      66.269  51.519  28.297  1.00 30.21           C  
ATOM    161  O3*   U A   8      65.321  51.887  27.314  1.00 32.41           O  
ATOM    162  C2*   U A   8      67.685  51.906  27.900  1.00 31.37           C  
ATOM    163  O2*   U A   8      67.743  53.224  27.433  1.00 27.02           O  
ATOM    164  C1*   U A   8      68.407  51.830  29.255  1.00 30.28           C  
ATOM    165  N1    U A   8      68.914  50.469  29.501  1.00 28.11           N  
ATOM    166  C2    U A   8      70.125  50.078  28.931  1.00 29.25           C  
ATOM    167  O2    U A   8      70.835  50.819  28.278  1.00 27.81           O  
ATOM    168  N3    U A   8      70.481  48.778  29.170  1.00 25.94           N  
ATOM    169  C4    U A   8      69.808  47.856  29.922  1.00 27.37           C  
ATOM    170  O4    U A   8      70.215  46.704  29.963  1.00 32.58           O  
ATOM    171  C5    U A   8      68.612  48.328  30.490  1.00 29.58           C  
ATOM    172  C6    U A   8      68.214  49.592  30.265  1.00 30.40           C  
ATOM    173  P     A A   9      64.755  50.731  26.311  1.00 32.05           P  
ATOM    174  O1P   A A   9      63.287  50.736  26.312  1.00 36.52           O  
ATOM    175  O2P   A A   9      65.447  49.415  26.519  1.00 31.33           O  
ATOM    176  O5*   A A   9      65.221  51.270  24.948  1.00 28.22           O  
ATOM    177  C5*   A A   9      64.646  52.433  24.369  1.00 32.91           C  
ATOM    178  C4*   A A   9      64.531  52.215  22.904  1.00 32.49           C  
ATOM    179  O4*   A A   9      65.887  52.090  22.406  1.00 35.09           O  
ATOM    180  C3*   A A   9      63.820  50.923  22.466  1.00 34.41           C  
ATOM    181  O3*   A A   9      63.140  51.180  21.236  1.00 36.11           O  
ATOM    182  C2*   A A   9      64.979  49.997  22.155  1.00 32.37           C  
ATOM    183  O2*   A A   9      64.686  49.016  21.194  1.00 35.87           O  
ATOM    184  C1*   A A   9      65.985  50.969  21.571  1.00 28.79           C  
ATOM    185  N9    A A   9      67.376  50.497  21.585  1.00 23.84           N  
ATOM    186  C8    A A   9      67.851  49.356  22.159  1.00 25.84           C  
ATOM    187  N7    A A   9      69.149  49.195  22.010  1.00 26.83           N  
ATOM    188  C5    A A   9      69.527  50.298  21.288  1.00 23.87           C  
ATOM    189  C6    A A   9      70.730  50.663  20.793  1.00 30.26           C  
ATOM    190  N6    A A   9      71.797  49.922  20.994  1.00 30.95           N  
ATOM    191  N1    A A   9      70.817  51.794  20.072  1.00 28.29           N  
ATOM    192  C2    A A   9      69.701  52.547  19.932  1.00 32.68           C  
ATOM    193  N3    A A   9      68.469  52.287  20.369  1.00 25.13           N  
ATOM    194  C4    A A   9      68.446  51.117  21.026  1.00 26.68           C  
HETATM  195  P   2MG A  10      61.504  51.328  21.232  1.00 44.21           P  
HETATM  196  O1P 2MG A  10      61.165  52.038  22.473  1.00 41.39           O  
HETATM  197  O2P 2MG A  10      61.216  51.946  19.892  1.00 41.97           O  
HETATM  198  O5* 2MG A  10      60.903  49.858  21.330  1.00 38.75           O  
HETATM  199  C5* 2MG A  10      59.437  49.660  21.397  1.00 42.74           C  
HETATM  200  C4* 2MG A  10      59.058  48.375  20.709  1.00 42.88           C  
HETATM  201  O4* 2MG A  10      59.575  48.416  19.351  1.00 44.02           O  
HETATM  202  C3* 2MG A  10      59.701  47.161  21.326  1.00 43.31           C  
HETATM  203  O3* 2MG A  10      58.874  46.647  22.357  1.00 45.12           O  
HETATM  204  C2* 2MG A  10      59.822  46.215  20.154  1.00 46.04           C  
HETATM  205  O2* 2MG A  10      58.533  45.637  19.943  1.00 47.96           O  
HETATM  206  C1* 2MG A  10      60.152  47.173  19.012  1.00 44.62           C  
HETATM  207  N9  2MG A  10      61.581  47.402  18.752  1.00 42.14           N  
HETATM  208  C8  2MG A  10      62.199  48.621  18.635  1.00 40.38           C  
HETATM  209  N7  2MG A  10      63.494  48.534  18.422  1.00 40.70           N  
HETATM  210  C5  2MG A  10      63.745  47.167  18.395  1.00 43.82           C  
HETATM  211  C6  2MG A  10      64.965  46.449  18.205  1.00 43.45           C  
HETATM  212  O6  2MG A  10      66.097  46.891  17.963  1.00 44.87           O  
HETATM  213  N1  2MG A  10      64.767  45.086  18.293  1.00 44.71           N  
HETATM  214  C2  2MG A  10      63.541  44.482  18.486  1.00 47.21           C  
HETATM  215  N2  2MG A  10      63.532  43.164  18.551  1.00 49.27           N  
HETATM  216  CM2 2MG A  10      62.220  42.454  18.591  1.00 52.10           C  
HETATM  217  N3  2MG A  10      62.411  45.125  18.614  1.00 45.85           N  
HETATM  218  C4  2MG A  10      62.574  46.451  18.582  1.00 43.27           C  
ATOM    219  P     C A  11      59.474  46.418  23.818  1.00 50.75           P  
ATOM    220  O1P   C A  11      58.367  46.417  24.802  1.00 49.46           O  
ATOM    221  O2P   C A  11      60.585  47.425  23.967  1.00 44.94           O  
ATOM    222  O5*   C A  11      60.064  44.937  23.797  1.00 49.65           O  
ATOM    223  C5*   C A  11      59.234  43.814  23.447  1.00 49.66           C  
ATOM    224  C4*   C A  11      60.091  42.608  23.221  1.00 50.13           C  
ATOM    225  O4*   C A  11      60.886  42.801  22.028  1.00 47.40           O  
ATOM    226  C3*   C A  11      61.091  42.406  24.335  1.00 52.64           C  
ATOM    227  O3*   C A  11      60.498  41.644  25.372  1.00 54.31           O  
ATOM    228  C2*   C A  11      62.252  41.701  23.640  1.00 51.51           C  
ATOM    229  O2*   C A  11      62.072  40.314  23.587  1.00 53.58           O  
ATOM    230  C1*   C A  11      62.189  42.294  22.230  1.00 48.88           C  
ATOM    231  N1    C A  11      63.145  43.397  21.999  1.00 46.21           N  
ATOM    232  C2    C A  11      64.484  43.091  21.738  1.00 45.74           C  
ATOM    233  O2    C A  11      64.833  41.895  21.708  1.00 47.49           O  
ATOM    234  N3    C A  11      65.365  44.106  21.527  1.00 42.07           N  
ATOM    235  C4    C A  11      64.941  45.376  21.555  1.00 42.69           C  
ATOM    236  N4    C A  11      65.829  46.353  21.301  1.00 38.00           N  
ATOM    237  C5    C A  11      63.586  45.709  21.822  1.00 43.75           C  
ATOM    238  C6    C A  11      62.732  44.698  22.035  1.00 45.15           C  
ATOM    239  P     U A  12      60.976  41.853  26.874  1.00 57.65           P  
ATOM    240  O1P   U A  12      60.123  41.039  27.783  1.00 59.26           O  
ATOM    241  O2P   U A  12      61.080  43.316  27.117  1.00 59.70           O  
ATOM    242  O5*   U A  12      62.441  41.244  26.886  1.00 55.93           O  
ATOM    243  C5*   U A  12      62.652  39.837  26.718  1.00 52.43           C  
ATOM    244  C4*   U A  12      64.121  39.544  26.594  1.00 48.54           C  
ATOM    245  O4*   U A  12      64.635  40.154  25.385  1.00 43.89           O  
ATOM    246  C3*   U A  12      65.015  40.119  27.684  1.00 50.30           C  
ATOM    247  O3*   U A  12      65.044  39.362  28.898  1.00 48.54           O  
ATOM    248  C2*   U A  12      66.384  40.159  27.007  1.00 48.05           C  
ATOM    249  O2*   U A  12      67.072  38.922  27.081  1.00 51.04           O  
ATOM    250  C1*   U A  12      66.013  40.478  25.565  1.00 47.40           C  
ATOM    251  N1    U A  12      66.260  41.889  25.203  1.00 42.96           N  
ATOM    252  C2    U A  12      67.565  42.250  24.872  1.00 43.65           C  
ATOM    253  O2    U A  12      68.512  41.462  24.914  1.00 41.62           O  
ATOM    254  N3    U A  12      67.736  43.577  24.520  1.00 45.49           N  
ATOM    255  C4    U A  12      66.755  44.567  24.497  1.00 41.58           C  
ATOM    256  O4    U A  12      67.069  45.722  24.203  1.00 41.28           O  
ATOM    257  C5    U A  12      65.436  44.114  24.864  1.00 40.07           C  
ATOM    258  C6    U A  12      65.241  42.831  25.193  1.00 44.32           C  
ATOM    259  P     C A  13      65.498  40.095  30.251  1.00 51.82           P  
ATOM    260  O1P   C A  13      65.477  39.159  31.399  1.00 52.57           O  
ATOM    261  O2P   C A  13      64.709  41.323  30.314  1.00 51.61           O  
ATOM    262  O5*   C A  13      67.010  40.399  29.846  1.00 48.38           O  
ATOM    263  C5*   C A  13      67.843  41.283  30.558  1.00 46.61           C  
ATOM    264  C4*   C A  13      69.102  41.516  29.762  1.00 43.48           C  
ATOM    265  O4*   C A  13      68.818  42.103  28.459  1.00 41.66           O  
ATOM    266  C3*   C A  13      70.084  42.477  30.389  1.00 46.01           C  
ATOM    267  O3*   C A  13      70.867  41.784  31.351  1.00 49.09           O  
ATOM    268  C2*   C A  13      70.925  42.922  29.204  1.00 42.77           C  
ATOM    269  O2*   C A  13      71.950  41.978  29.004  1.00 45.11           O  
ATOM    270  C1*   C A  13      69.909  42.919  28.054  1.00 39.03           C  
ATOM    271  N1    C A  13      69.390  44.256  27.668  1.00 36.14           N  
ATOM    272  C2    C A  13      70.236  45.162  26.999  1.00 31.39           C  
ATOM    273  O2    C A  13      71.429  44.834  26.781  1.00 33.97           O  
ATOM    274  N3    C A  13      69.736  46.376  26.604  1.00 30.71           N  
ATOM    275  C4    C A  13      68.478  46.702  26.907  1.00 27.53           C  
ATOM    276  N4    C A  13      68.050  47.913  26.576  1.00 29.33           N  
ATOM    277  C5    C A  13      67.596  45.795  27.581  1.00 30.21           C  
ATOM    278  C6    C A  13      68.085  44.602  27.937  1.00 31.74           C  
ATOM    279  P     A A  14      71.499  42.582  32.585  1.00 52.69           P  
ATOM    280  O1P   A A  14      71.592  41.647  33.723  1.00 56.81           O  
ATOM    281  O2P   A A  14      70.795  43.877  32.732  1.00 54.35           O  
ATOM    282  O5*   A A  14      72.996  42.894  32.143  1.00 52.44           O  
ATOM    283  C5*   A A  14      73.291  44.004  31.337  1.00 44.62           C  
ATOM    284  C4*   A A  14      74.612  43.815  30.626  1.00 38.21           C  
ATOM    285  O4*   A A  14      74.372  44.178  29.229  1.00 37.15           O  
ATOM    286  C3*   A A  14      75.617  44.841  31.120  1.00 39.31           C  
ATOM    287  O3*   A A  14      76.409  44.373  32.214  1.00 35.14           O  
ATOM    288  C2*   A A  14      76.410  45.187  29.878  1.00 35.42           C  
ATOM    289  O2*   A A  14      77.406  44.222  29.562  1.00 37.00           O  
ATOM    290  C1*   A A  14      75.325  45.147  28.805  1.00 32.12           C  
ATOM    291  N9    A A  14      74.639  46.437  28.568  1.00 30.32           N  
ATOM    292  C8    A A  14      73.332  46.800  28.850  1.00 28.83           C  
ATOM    293  N7    A A  14      73.030  48.029  28.495  1.00 27.89           N  
ATOM    294  C5    A A  14      74.205  48.496  27.963  1.00 28.40           C  
ATOM    295  C6    A A  14      74.551  49.722  27.451  1.00 28.49           C  
ATOM    296  N6    A A  14      73.715  50.778  27.422  1.00 27.62           N  
ATOM    297  N1    A A  14      75.820  49.867  26.972  1.00 28.90           N  
ATOM    298  C2    A A  14      76.658  48.824  27.058  1.00 26.09           C  
ATOM    299  N3    A A  14      76.449  47.633  27.546  1.00 34.21           N  
ATOM    300  C4    A A  14      75.194  47.523  27.993  1.00 28.42           C  
ATOM    301  P     G A  15      76.463  45.227  33.560  1.00 38.35           P  
ATOM    302  O1P   G A  15      77.577  44.561  34.373  1.00 36.34           O  
ATOM    303  O2P   G A  15      75.020  45.308  34.127  1.00 36.07           O  
ATOM    304  O5*   G A  15      76.977  46.682  33.179  1.00 34.17           O  
ATOM    305  C5*   G A  15      78.216  46.873  32.475  1.00 37.95           C  
ATOM    306  C4*   G A  15      78.274  48.248  31.867  1.00 34.05           C  
ATOM    307  O4*   G A  15      77.353  48.424  30.762  1.00 36.42           O  
ATOM    308  C3*   G A  15      77.992  49.400  32.793  1.00 41.25           C  
ATOM    309  O3*   G A  15      79.176  49.641  33.526  1.00 55.04           O  
ATOM    310  C2*   G A  15      77.696  50.528  31.803  1.00 37.22           C  
ATOM    311  O2*   G A  15      78.880  51.102  31.276  1.00 33.49           O  
ATOM    312  C1*   G A  15      76.941  49.779  30.686  1.00 33.07           C  
ATOM    313  N9    G A  15      75.505  49.813  30.961  1.00 30.40           N  
ATOM    314  C8    G A  15      74.775  48.836  31.612  1.00 30.39           C  
ATOM    315  N7    G A  15      73.537  49.203  31.861  1.00 28.29           N  
ATOM    316  C5    G A  15      73.439  50.464  31.298  1.00 27.38           C  
ATOM    317  C6    G A  15      72.351  51.383  31.271  1.00 24.88           C  
ATOM    318  O6    G A  15      71.261  51.260  31.756  1.00 26.54           O  
ATOM    319  N1    G A  15      72.683  52.569  30.607  1.00 23.73           N  
ATOM    320  C2    G A  15      73.896  52.858  30.094  1.00 29.75           C  
ATOM    321  N2    G A  15      74.047  54.083  29.581  1.00 27.02           N  
ATOM    322  N3    G A  15      74.925  52.008  30.089  1.00 28.46           N  
ATOM    323  C4    G A  15      74.632  50.842  30.714  1.00 29.28           C  
HETATM  324  P   H2U A  16      79.106  49.914  35.099  1.00 64.01           P  
HETATM  325  O1P H2U A  16      77.803  50.520  35.400  1.00 58.28           O  
HETATM  326  O2P H2U A  16      79.533  48.676  35.816  1.00 67.91           O  
HETATM  327  O5* H2U A  16      80.270  50.994  35.265  1.00 70.49           O  
HETATM  328  C5* H2U A  16      81.110  51.317  34.115  1.00 77.82           C  
HETATM  329  C4* H2U A  16      80.514  52.486  33.353  1.00 82.34           C  
HETATM  330  O4* H2U A  16      79.081  52.313  33.356  1.00 85.70           O  
HETATM  331  C3* H2U A  16      80.758  53.821  34.030  1.00 84.30           C  
HETATM  332  O3* H2U A  16      81.907  54.422  33.414  1.00 84.12           O  
HETATM  333  C1* H2U A  16      78.428  53.548  33.551  1.00 88.13           C  
HETATM  334  C2* H2U A  16      79.505  54.639  33.690  1.00 86.71           C  
HETATM  335  O2* H2U A  16      79.637  55.391  32.493  1.00 88.25           O  
HETATM  336  N1  H2U A  16      77.347  53.323  34.582  1.00 91.19           N  
HETATM  337  C2  H2U A  16      76.119  52.865  34.160  1.00 92.39           C  
HETATM  338  O2  H2U A  16      75.885  52.463  33.033  1.00 92.20           O  
HETATM  339  N3  H2U A  16      75.123  52.894  35.107  1.00 93.28           N  
HETATM  340  C4  H2U A  16      75.289  52.711  36.458  1.00 93.34           C  
HETATM  341  O4  H2U A  16      74.309  52.695  37.208  1.00 92.66           O  
HETATM  342  C5  H2U A  16      76.696  52.479  36.909  1.00 93.77           C  
HETATM  343  C6  H2U A  16      77.717  53.238  36.039  1.00 93.22           C  
HETATM  344  P   H2U A  17      83.371  53.708  33.472  1.00 82.84           P  
HETATM  345  O1P H2U A  17      83.746  53.377  32.068  1.00 83.70           O  
HETATM  346  O2P H2U A  17      83.498  52.655  34.529  1.00 83.82           O  
HETATM  347  O5* H2U A  17      84.277  54.923  33.943  1.00 81.72           O  
HETATM  348  C5* H2U A  17      83.692  55.978  34.736  1.00 76.14           C  
HETATM  349  C4* H2U A  17      84.176  55.886  36.150  1.00 71.56           C  
HETATM  350  O4* H2U A  17      85.622  55.872  36.137  1.00 71.61           O  
HETATM  351  C3* H2U A  17      83.738  57.031  37.055  1.00 67.99           C  
HETATM  352  O3* H2U A  17      82.553  56.582  37.718  1.00 60.02           O  
HETATM  353  C1* H2U A  17      86.102  56.903  36.958  1.00 71.64           C  
HETATM  354  C2* H2U A  17      84.964  57.213  37.948  1.00 71.27           C  
HETATM  355  O2* H2U A  17      85.004  56.273  39.021  1.00 73.23           O  
HETATM  356  N1  H2U A  17      86.579  57.954  36.004  1.00 72.27           N  
HETATM  357  C2  H2U A  17      87.702  58.662  36.301  1.00 71.21           C  
HETATM  358  O2  H2U A  17      87.834  59.359  37.287  1.00 72.68           O  
HETATM  359  N3  H2U A  17      88.693  58.585  35.358  1.00 69.04           N  
HETATM  360  C4  H2U A  17      88.711  57.779  34.244  1.00 68.89           C  
HETATM  361  O4  H2U A  17      89.766  57.616  33.620  1.00 64.81           O  
HETATM  362  C5  H2U A  17      87.401  57.154  33.864  1.00 69.36           C  
HETATM  363  C6  H2U A  17      86.257  57.828  34.577  1.00 71.72           C  
ATOM    364  P     G A  18      81.804  57.491  38.803  1.00 53.35           P  
ATOM    365  O1P   G A  18      82.773  58.100  39.715  1.00 56.70           O  
ATOM    366  O2P   G A  18      80.724  56.638  39.368  1.00 56.95           O  
ATOM    367  O5*   G A  18      81.038  58.580  37.950  1.00 45.07           O  
ATOM    368  C5*   G A  18      80.288  58.201  36.778  1.00 37.17           C  
ATOM    369  C4*   G A  18      80.100  59.412  35.902  1.00 33.24           C  
ATOM    370  O4*   G A  18      79.417  60.430  36.705  1.00 29.49           O  
ATOM    371  C3*   G A  18      81.426  60.038  35.456  1.00 29.19           C  
ATOM    372  O3*   G A  18      81.313  60.691  34.173  1.00 27.98           O  
ATOM    373  C2*   G A  18      81.638  61.165  36.437  1.00 26.61           C  
ATOM    374  O2*   G A  18      82.417  62.205  35.773  1.00 31.88           O  
ATOM    375  C1*   G A  18      80.191  61.615  36.658  1.00 31.24           C  
ATOM    376  N9    G A  18      79.893  62.457  37.818  1.00 25.73           N  
ATOM    377  C8    G A  18      80.399  62.361  39.094  1.00 29.39           C  
ATOM    378  N7    G A  18      79.992  63.332  39.883  1.00 29.66           N  
ATOM    379  C5    G A  18      79.165  64.095  39.074  1.00 26.03           C  
ATOM    380  C6    G A  18      78.469  65.259  39.359  1.00 29.73           C  
ATOM    381  O6    G A  18      78.491  65.889  40.411  1.00 31.10           O  
ATOM    382  N1    G A  18      77.698  65.711  38.257  1.00 25.44           N  
ATOM    383  C2    G A  18      77.634  65.077  37.076  1.00 25.15           C  
ATOM    384  N2    G A  18      76.850  65.605  36.127  1.00 28.06           N  
ATOM    385  N3    G A  18      78.312  63.960  36.797  1.00 31.81           N  
ATOM    386  C4    G A  18      79.053  63.539  37.817  1.00 27.46           C  
ATOM    387  P     G A  19      81.705  59.935  32.855  1.00 34.61           P  
ATOM    388  O1P   G A  19      80.751  58.780  32.632  1.00 31.49           O  
ATOM    389  O2P   G A  19      83.185  59.683  32.792  1.00 35.14           O  
ATOM    390  O5*   G A  19      81.429  61.065  31.780  1.00 30.92           O  
ATOM    391  C5*   G A  19      80.053  61.459  31.456  1.00 35.49           C  
ATOM    392  C4*   G A  19      80.105  62.508  30.407  1.00 32.28           C  
ATOM    393  O4*   G A  19      80.779  63.631  30.991  1.00 33.92           O  
ATOM    394  C3*   G A  19      80.907  62.116  29.171  1.00 34.56           C  
ATOM    395  O3*   G A  19      80.389  62.868  28.083  1.00 33.76           O  
ATOM    396  C2*   G A  19      82.305  62.679  29.462  1.00 35.83           C  
ATOM    397  O2*   G A  19      82.892  63.160  28.284  1.00 36.38           O  
ATOM    398  C1*   G A  19      81.965  63.912  30.283  1.00 31.26           C  
ATOM    399  N9    G A  19      82.922  64.301  31.303  1.00 34.09           N  
ATOM    400  C8    G A  19      83.808  63.510  32.002  1.00 34.35           C  
ATOM    401  N7    G A  19      84.330  64.135  33.026  1.00 29.97           N  
ATOM    402  C5    G A  19      83.803  65.412  32.954  1.00 30.15           C  
ATOM    403  C6    G A  19      83.998  66.553  33.810  1.00 30.28           C  
ATOM    404  O6    G A  19      84.715  66.645  34.861  1.00 30.49           O  
ATOM    405  N1    G A  19      83.246  67.654  33.362  1.00 30.06           N  
ATOM    406  C2    G A  19      82.438  67.667  32.262  1.00 28.68           C  
ATOM    407  N2    G A  19      81.753  68.812  32.044  1.00 31.69           N  
ATOM    408  N3    G A  19      82.279  66.621  31.443  1.00 33.35           N  
ATOM    409  C4    G A  19      82.976  65.537  31.862  1.00 29.77           C  
ATOM    410  P     G A  20      79.212  62.270  27.169  1.00 37.57           P  
ATOM    411  O1P   G A  20      78.632  63.435  26.460  1.00 37.48           O  
ATOM    412  O2P   G A  20      78.418  61.535  28.129  1.00 34.93           O  
ATOM    413  O5*   G A  20      79.952  61.358  26.110  1.00 35.27           O  
ATOM    414  C5*   G A  20      80.561  61.904  24.945  1.00 34.78           C  
ATOM    415  C4*   G A  20      80.529  60.879  23.849  1.00 33.56           C  
ATOM    416  O4*   G A  20      81.420  59.808  24.203  1.00 34.69           O  
ATOM    417  C3*   G A  20      79.164  60.229  23.679  1.00 34.55           C  
ATOM    418  O3*   G A  20      78.373  60.964  22.748  1.00 32.07           O  
ATOM    419  C2*   G A  20      79.540  58.879  23.115  1.00 31.94           C  
ATOM    420  O2*   G A  20      79.935  59.005  21.776  1.00 29.14           O  
ATOM    421  C1*   G A  20      80.836  58.580  23.876  1.00 33.29           C  
ATOM    422  N9    G A  20      80.557  57.879  25.128  1.00 28.16           N  
ATOM    423  C8    G A  20      80.808  58.309  26.413  1.00 28.17           C  
ATOM    424  N7    G A  20      80.362  57.471  27.316  1.00 27.62           N  
ATOM    425  C5    G A  20      79.819  56.431  26.576  1.00 26.85           C  
ATOM    426  C6    G A  20      79.214  55.251  26.985  1.00 29.07           C  
ATOM    427  O6    G A  20      79.010  54.845  28.141  1.00 30.12           O  
ATOM    428  N1    G A  20      78.811  54.474  25.894  1.00 26.18           N  
ATOM    429  C2    G A  20      78.993  54.795  24.576  1.00 28.16           C  
ATOM    430  N2    G A  20      78.546  53.879  23.648  1.00 24.71           N  
ATOM    431  N3    G A  20      79.567  55.902  24.177  1.00 25.02           N  
ATOM    432  C4    G A  20      79.951  56.671  25.225  1.00 27.86           C  
ATOM    433  P     A A  21      76.960  61.541  23.213  1.00 35.97           P  
ATOM    434  O1P   A A  21      76.324  62.054  21.989  1.00 35.07           O  
ATOM    435  O2P   A A  21      77.193  62.456  24.350  1.00 35.45           O  
ATOM    436  O5*   A A  21      76.166  60.289  23.790  1.00 34.57           O  
ATOM    437  C5*   A A  21      75.604  59.329  22.914  1.00 33.17           C  
ATOM    438  C4*   A A  21      75.622  58.001  23.580  1.00 34.52           C  
ATOM    439  O4*   A A  21      74.864  57.958  24.808  1.00 29.25           O  
ATOM    440  C3*   A A  21      75.125  56.852  22.735  1.00 33.11           C  
ATOM    441  O3*   A A  21      76.250  56.581  21.883  1.00 36.69           O  
ATOM    442  C2*   A A  21      74.815  55.806  23.788  1.00 33.65           C  
ATOM    443  O2*   A A  21      76.034  55.158  24.220  1.00 30.39           O  
ATOM    444  C1*   A A  21      74.304  56.666  24.955  1.00 30.28           C  
ATOM    445  N9    A A  21      72.843  56.834  25.065  1.00 28.13           N  
ATOM    446  C8    A A  21      72.122  57.968  24.720  1.00 24.90           C  
ATOM    447  N7    A A  21      70.828  57.864  24.968  1.00 27.77           N  
ATOM    448  C5    A A  21      70.695  56.581  25.516  1.00 25.63           C  
ATOM    449  C6    A A  21      69.596  55.872  26.031  1.00 28.10           C  
ATOM    450  N6    A A  21      68.327  56.341  26.073  1.00 24.98           N  
ATOM    451  N1    A A  21      69.817  54.668  26.536  1.00 27.40           N  
ATOM    452  C2    A A  21      71.064  54.186  26.531  1.00 26.83           C  
ATOM    453  N3    A A  21      72.169  54.734  26.097  1.00 26.24           N  
ATOM    454  C4    A A  21      71.927  55.948  25.580  1.00 28.94           C  
ATOM    455  P     G A  22      76.122  55.607  20.624  1.00 36.77           P  
ATOM    456  O1P   G A  22      77.347  55.811  19.853  1.00 34.15           O  
ATOM    457  O2P   G A  22      74.796  55.966  20.020  1.00 38.21           O  
ATOM    458  O5*   G A  22      76.107  54.255  21.420  1.00 39.67           O  
ATOM    459  C5*   G A  22      75.588  53.058  20.896  1.00 33.77           C  
ATOM    460  C4*   G A  22      76.292  51.916  21.581  1.00 33.57           C  
ATOM    461  O4*   G A  22      76.032  51.934  23.018  1.00 28.98           O  
ATOM    462  C3*   G A  22      75.794  50.585  21.102  1.00 31.27           C  
ATOM    463  O3*   G A  22      76.427  50.216  19.874  1.00 36.07           O  
ATOM    464  C2*   G A  22      75.986  49.707  22.315  1.00 30.09           C  
ATOM    465  O2*   G A  22      77.321  49.252  22.478  1.00 29.86           O  
ATOM    466  C1*   G A  22      75.605  50.671  23.444  1.00 25.75           C  
ATOM    467  N9    G A  22      74.157  50.722  23.757  1.00 27.25           N  
ATOM    468  C8    G A  22      73.306  51.785  23.618  1.00 27.24           C  
ATOM    469  N7    G A  22      72.095  51.529  24.052  1.00 28.94           N  
ATOM    470  C5    G A  22      72.147  50.220  24.476  1.00 28.23           C  
ATOM    471  C6    G A  22      71.138  49.385  25.011  1.00 28.84           C  
ATOM    472  O6    G A  22      69.936  49.648  25.252  1.00 29.10           O  
ATOM    473  N1    G A  22      71.624  48.128  25.270  1.00 27.83           N  
ATOM    474  C2    G A  22      72.899  47.703  25.050  1.00 30.21           C  
ATOM    475  N2    G A  22      73.158  46.409  25.364  1.00 27.00           N  
ATOM    476  N3    G A  22      73.837  48.454  24.566  1.00 29.25           N  
ATOM    477  C4    G A  22      73.400  49.695  24.288  1.00 25.81           C  
ATOM    478  P     A A  23      75.571  49.404  18.784  1.00 38.64           P  
ATOM    479  O1P   A A  23      76.230  49.559  17.497  1.00 35.37           O  
ATOM    480  O2P   A A  23      74.080  49.759  18.881  1.00 35.35           O  
ATOM    481  O5*   A A  23      75.636  47.938  19.332  1.00 37.83           O  
ATOM    482  C5*   A A  23      76.912  47.318  19.633  1.00 37.02           C  
ATOM    483  C4*   A A  23      76.705  45.975  20.257  1.00 39.13           C  
ATOM    484  O4*   A A  23      76.089  46.103  21.547  1.00 38.15           O  
ATOM    485  C3*   A A  23      75.794  45.028  19.486  1.00 40.14           C  
ATOM    486  O3*   A A  23      76.500  44.368  18.446  1.00 42.27           O  
ATOM    487  C2*   A A  23      75.356  44.060  20.563  1.00 39.45           C  
ATOM    488  O2*   A A  23      76.423  43.125  20.754  1.00 42.00           O  
ATOM    489  C1*   A A  23      75.210  45.010  21.771  1.00 37.98           C  
ATOM    490  N9    A A  23      73.858  45.578  21.885  1.00 35.38           N  
ATOM    491  C8    A A  23      73.461  46.822  21.439  1.00 35.25           C  
ATOM    492  N7    A A  23      72.234  47.140  21.772  1.00 33.59           N  
ATOM    493  C5    A A  23      71.772  46.017  22.438  1.00 34.26           C  
ATOM    494  C6    A A  23      70.529  45.710  23.019  1.00 32.33           C  
ATOM    495  N6    A A  23      69.521  46.594  23.138  1.00 34.94           N  
ATOM    496  N1    A A  23      70.368  44.470  23.517  1.00 35.35           N  
ATOM    497  C2    A A  23      71.405  43.615  23.470  1.00 31.66           C  
ATOM    498  N3    A A  23      72.643  43.811  23.008  1.00 32.48           N  
ATOM    499  C4    A A  23      72.758  45.035  22.487  1.00 34.85           C  
ATOM    500  P     G A  24      75.691  43.814  17.170  1.00 44.02           P  
ATOM    501  O1P   G A  24      76.744  43.269  16.286  1.00 46.42           O  
ATOM    502  O2P   G A  24      74.752  44.889  16.691  1.00 43.05           O  
ATOM    503  O5*   G A  24      74.795  42.670  17.756  1.00 40.73           O  
ATOM    504  C5*   G A  24      75.378  41.492  18.218  1.00 46.30           C  
ATOM    505  C4*   G A  24      74.313  40.622  18.747  1.00 47.22           C  
ATOM    506  O4*   G A  24      73.799  41.198  19.975  1.00 46.23           O  
ATOM    507  C3*   G A  24      73.094  40.484  17.855  1.00 47.23           C  
ATOM    508  O3*   G A  24      73.287  39.486  16.850  1.00 51.92           O  
ATOM    509  C2*   G A  24      72.056  40.037  18.867  1.00 47.49           C  
ATOM    510  O2*   G A  24      72.324  38.676  19.169  1.00 46.29           O  
ATOM    511  C1*   G A  24      72.412  40.925  20.073  1.00 45.11           C  
ATOM    512  N9    G A  24      71.687  42.195  20.013  1.00 42.04           N  
ATOM    513  C8    G A  24      72.126  43.377  19.471  1.00 42.70           C  
ATOM    514  N7    G A  24      71.218  44.315  19.472  1.00 40.83           N  
ATOM    515  C5    G A  24      70.114  43.715  20.070  1.00 42.41           C  
ATOM    516  C6    G A  24      68.831  44.242  20.336  1.00 44.36           C  
ATOM    517  O6    G A  24      68.390  45.396  20.057  1.00 44.29           O  
ATOM    518  N1    G A  24      68.011  43.299  20.965  1.00 42.61           N  
ATOM    519  C2    G A  24      68.390  42.013  21.284  1.00 44.44           C  
ATOM    520  N2    G A  24      67.446  41.234  21.843  1.00 42.33           N  
ATOM    521  N3    G A  24      69.599  41.521  21.050  1.00 41.16           N  
ATOM    522  C4    G A  24      70.396  42.419  20.435  1.00 42.35           C  
ATOM    523  P     C A  25      72.660  39.699  15.392  1.00 49.58           P  
ATOM    524  O1P   C A  25      73.082  38.530  14.539  1.00 54.65           O  
ATOM    525  O2P   C A  25      72.897  41.067  14.886  1.00 47.78           O  
ATOM    526  O5*   C A  25      71.110  39.469  15.645  1.00 48.69           O  
ATOM    527  C5*   C A  25      70.619  38.180  16.052  1.00 50.80           C  
ATOM    528  C4*   C A  25      69.170  38.275  16.491  1.00 50.78           C  
ATOM    529  O4*   C A  25      69.098  39.152  17.653  1.00 55.06           O  
ATOM    530  C3*   C A  25      68.156  38.899  15.531  1.00 53.87           C  
ATOM    531  O3*   C A  25      67.629  38.012  14.533  1.00 52.62           O  
ATOM    532  C2*   C A  25      67.045  39.336  16.481  1.00 54.24           C  
ATOM    533  O2*   C A  25      66.276  38.238  16.922  1.00 56.92           O  
ATOM    534  C1*   C A  25      67.844  39.818  17.685  1.00 51.93           C  
ATOM    535  N1    C A  25      68.040  41.263  17.589  1.00 51.02           N  
ATOM    536  C2    C A  25      67.005  42.092  18.016  1.00 51.02           C  
ATOM    537  O2    C A  25      65.989  41.579  18.475  1.00 49.40           O  
ATOM    538  N3    C A  25      67.136  43.428  17.897  1.00 50.42           N  
ATOM    539  C4    C A  25      68.240  43.949  17.341  1.00 50.86           C  
ATOM    540  N4    C A  25      68.316  45.287  17.212  1.00 46.56           N  
ATOM    541  C5    C A  25      69.316  43.129  16.899  1.00 48.98           C  
ATOM    542  C6    C A  25      69.185  41.802  17.062  1.00 50.79           C  
HETATM  543  P   M2G A  26      67.297  38.593  13.077  1.00 53.21           P  
HETATM  544  O1P M2G A  26      68.469  39.430  12.691  1.00 53.23           O  
HETATM  545  O2P M2G A  26      66.835  37.487  12.148  1.00 56.89           O  
HETATM  546  O5* M2G A  26      66.058  39.572  13.258  1.00 49.44           O  
HETATM  547  C5* M2G A  26      64.865  39.112  13.893  1.00 49.64           C  
HETATM  548  C4* M2G A  26      63.938  40.267  14.143  1.00 49.36           C  
HETATM  549  O4* M2G A  26      64.443  41.119  15.209  1.00 48.17           O  
HETATM  550  C3* M2G A  26      63.719  41.196  12.968  1.00 48.54           C  
HETATM  551  O3* M2G A  26      62.681  40.699  12.152  1.00 52.14           O  
HETATM  552  C2* M2G A  26      63.273  42.477  13.644  1.00 47.40           C  
HETATM  553  O2* M2G A  26      61.905  42.387  14.025  1.00 47.93           O  
HETATM  554  C1* M2G A  26      64.147  42.470  14.901  1.00 46.92           C  
HETATM  555  N9  M2G A  26      65.403  43.195  14.677  1.00 41.61           N  
HETATM  556  C8  M2G A  26      66.657  42.687  14.400  1.00 45.40           C  
HETATM  557  N7  M2G A  26      67.543  43.618  14.188  1.00 44.60           N  
HETATM  558  C5  M2G A  26      66.836  44.807  14.371  1.00 43.23           C  
HETATM  559  C6  M2G A  26      67.253  46.174  14.285  1.00 41.77           C  
HETATM  560  O6  M2G A  26      68.372  46.637  14.008  1.00 45.23           O  
HETATM  561  N1  M2G A  26      66.209  47.048  14.544  1.00 42.85           N  
HETATM  562  C2  M2G A  26      64.926  46.678  14.840  1.00 42.16           C  
HETATM  563  N2  M2G A  26      64.015  47.663  15.061  1.00 41.43           N  
HETATM  564  N3  M2G A  26      64.531  45.410  14.927  1.00 42.18           N  
HETATM  565  C4  M2G A  26      65.524  44.546  14.680  1.00 43.64           C  
HETATM  566  CM1 M2G A  26      64.404  49.075  15.158  1.00 44.51           C  
HETATM  567  CM2 M2G A  26      62.594  47.288  15.283  1.00 41.28           C  
ATOM    568  P     C A  27      62.860  40.709  10.569  1.00 53.96           P  
ATOM    569  O1P   C A  27      61.690  39.909  10.056  1.00 58.78           O  
ATOM    570  O2P   C A  27      64.244  40.343  10.165  1.00 54.14           O  
ATOM    571  O5*   C A  27      62.614  42.218  10.155  1.00 53.83           O  
ATOM    572  C5*   C A  27      61.437  42.899  10.561  1.00 51.71           C  
ATOM    573  C4*   C A  27      61.594  44.378  10.321  1.00 53.47           C  
ATOM    574  O4*   C A  27      62.476  44.955  11.330  1.00 51.31           O  
ATOM    575  C3*   C A  27      62.210  44.814   9.004  1.00 54.38           C  
ATOM    576  O3*   C A  27      61.256  44.865   7.943  1.00 60.47           O  
ATOM    577  C2*   C A  27      62.693  46.223   9.339  1.00 51.78           C  
ATOM    578  O2*   C A  27      61.640  47.167   9.301  1.00 51.18           O  
ATOM    579  C1*   C A  27      63.159  46.064  10.787  1.00 48.68           C  
ATOM    580  N1    C A  27      64.612  45.818  10.877  1.00 42.59           N  
ATOM    581  C2    C A  27      65.472  46.868  10.634  1.00 44.48           C  
ATOM    582  O2    C A  27      64.981  47.978  10.348  1.00 42.73           O  
ATOM    583  N3    C A  27      66.821  46.659  10.722  1.00 42.28           N  
ATOM    584  C4    C A  27      67.275  45.452  11.056  1.00 43.75           C  
ATOM    585  N4    C A  27      68.586  45.272  11.180  1.00 44.57           N  
ATOM    586  C5    C A  27      66.402  44.364  11.291  1.00 44.20           C  
ATOM    587  C6    C A  27      65.095  44.589  11.192  1.00 44.33           C  
ATOM    588  P     C A  28      61.715  44.551   6.429  1.00 61.60           P  
ATOM    589  O1P   C A  28      60.464  44.441   5.640  1.00 62.15           O  
ATOM    590  O2P   C A  28      62.648  43.401   6.479  1.00 61.61           O  
ATOM    591  O5*   C A  28      62.553  45.816   5.941  1.00 56.70           O  
ATOM    592  C5*   C A  28      61.910  47.047   5.647  1.00 57.75           C  
ATOM    593  C4*   C A  28      62.918  48.099   5.246  1.00 58.29           C  
ATOM    594  O4*   C A  28      63.719  48.476   6.402  1.00 56.80           O  
ATOM    595  C3*   C A  28      63.942  47.721   4.187  1.00 60.32           C  
ATOM    596  O3*   C A  28      63.466  47.779   2.838  1.00 64.19           O  
ATOM    597  C2*   C A  28      65.056  48.727   4.461  1.00 57.91           C  
ATOM    598  O2*   C A  28      64.774  50.025   3.952  1.00 56.60           O  
ATOM    599  C1*   C A  28      65.046  48.767   5.987  1.00 56.20           C  
ATOM    600  N1    C A  28      65.958  47.725   6.532  1.00 51.41           N  
ATOM    601  C2    C A  28      67.338  47.999   6.616  1.00 49.26           C  
ATOM    602  O2    C A  28      67.761  49.119   6.240  1.00 45.60           O  
ATOM    603  N3    C A  28      68.175  47.034   7.101  1.00 50.25           N  
ATOM    604  C4    C A  28      67.676  45.854   7.493  1.00 50.21           C  
ATOM    605  N4    C A  28      68.515  44.934   7.978  1.00 50.49           N  
ATOM    606  C5    C A  28      66.284  45.560   7.412  1.00 50.83           C  
ATOM    607  C6    C A  28      65.473  46.516   6.936  1.00 50.67           C  
ATOM    608  P     A A  29      64.134  46.810   1.715  1.00 65.77           P  
ATOM    609  O1P   A A  29      63.389  47.109   0.457  1.00 67.04           O  
ATOM    610  O2P   A A  29      64.179  45.395   2.219  1.00 63.95           O  
ATOM    611  O5*   A A  29      65.641  47.328   1.540  1.00 62.83           O  
ATOM    612  C5*   A A  29      65.903  48.647   0.972  1.00 63.87           C  
ATOM    613  C4*   A A  29      67.384  49.009   1.027  1.00 63.27           C  
ATOM    614  O4*   A A  29      67.815  49.121   2.412  1.00 62.56           O  
ATOM    615  C3*   A A  29      68.423  48.084   0.399  1.00 63.00           C  
ATOM    616  O3*   A A  29      68.576  48.226  -1.003  1.00 65.01           O  
ATOM    617  C2*   A A  29      69.702  48.526   1.092  1.00 61.38           C  
ATOM    618  O2*   A A  29      70.234  49.721   0.530  1.00 61.02           O  
ATOM    619  C1*   A A  29      69.196  48.789   2.509  1.00 59.63           C  
ATOM    620  N9    A A  29      69.335  47.576   3.315  1.00 54.19           N  
ATOM    621  C8    A A  29      68.378  46.647   3.653  1.00 54.48           C  
ATOM    622  N7    A A  29      68.847  45.630   4.321  1.00 53.76           N  
ATOM    623  C5    A A  29      70.198  45.923   4.448  1.00 52.14           C  
ATOM    624  C6    A A  29      71.245  45.231   5.038  1.00 50.79           C  
ATOM    625  N6    A A  29      71.082  44.066   5.687  1.00 50.54           N  
ATOM    626  N1    A A  29      72.481  45.776   4.952  1.00 48.83           N  
ATOM    627  C2    A A  29      72.620  46.942   4.329  1.00 48.49           C  
ATOM    628  N3    A A  29      71.699  47.698   3.751  1.00 51.37           N  
ATOM    629  C4    A A  29      70.502  47.119   3.845  1.00 50.75           C  
ATOM    630  P     G A  30      69.148  46.986  -1.875  1.00 63.15           P  
ATOM    631  O1P   G A  30      68.830  47.365  -3.277  1.00 64.99           O  
ATOM    632  O2P   G A  30      68.659  45.688  -1.323  1.00 61.44           O  
ATOM    633  O5*   G A  30      70.738  46.979  -1.648  1.00 64.00           O  
ATOM    634  C5*   G A  30      71.552  48.159  -1.889  1.00 63.25           C  
ATOM    635  C4*   G A  30      72.933  47.996  -1.281  1.00 63.56           C  
ATOM    636  O4*   G A  30      72.830  47.805   0.157  1.00 62.05           O  
ATOM    637  C3*   G A  30      73.716  46.785  -1.754  1.00 65.21           C  
ATOM    638  O3*   G A  30      74.380  47.055  -2.978  1.00 68.93           O  
ATOM    639  C2*   G A  30      74.690  46.532  -0.603  1.00 63.90           C  
ATOM    640  O2*   G A  30      75.815  47.379  -0.636  1.00 66.48           O  
ATOM    641  C1*   G A  30      73.839  46.903   0.605  1.00 60.48           C  
ATOM    642  N9    G A  30      73.204  45.719   1.184  1.00 56.62           N  
ATOM    643  C8    G A  30      71.872  45.376   1.155  1.00 55.59           C  
ATOM    644  N7    G A  30      71.623  44.235   1.750  1.00 54.28           N  
ATOM    645  C5    G A  30      72.868  43.798   2.203  1.00 52.26           C  
ATOM    646  C6    G A  30      73.238  42.611   2.909  1.00 51.92           C  
ATOM    647  O6    G A  30      72.512  41.677   3.315  1.00 48.72           O  
ATOM    648  N1    G A  30      74.614  42.562   3.133  1.00 53.22           N  
ATOM    649  C2    G A  30      75.513  43.524   2.738  1.00 52.55           C  
ATOM    650  N2    G A  30      76.797  43.282   3.015  1.00 52.61           N  
ATOM    651  N3    G A  30      75.180  44.635   2.106  1.00 52.73           N  
ATOM    652  C4    G A  30      73.853  44.703   1.866  1.00 53.73           C  
ATOM    653  P     A A  31      74.530  45.886  -4.075  1.00 69.42           P  
ATOM    654  O1P   A A  31      75.319  46.471  -5.195  1.00 71.00           O  
ATOM    655  O2P   A A  31      73.193  45.303  -4.337  1.00 68.84           O  
ATOM    656  O5*   A A  31      75.405  44.785  -3.333  1.00 67.81           O  
ATOM    657  C5*   A A  31      76.782  45.018  -3.047  1.00 66.34           C  
ATOM    658  C4*   A A  31      77.346  43.848  -2.300  1.00 65.86           C  
ATOM    659  O4*   A A  31      76.688  43.760  -1.020  1.00 64.48           O  
ATOM    660  C3*   A A  31      77.123  42.485  -2.938  1.00 67.80           C  
ATOM    661  O3*   A A  31      78.102  42.185  -3.933  1.00 69.43           O  
ATOM    662  C2*   A A  31      77.166  41.543  -1.735  1.00 66.49           C  
ATOM    663  O2*   A A  31      78.476  41.202  -1.324  1.00 70.09           O  
ATOM    664  C1*   A A  31      76.518  42.403  -0.649  1.00 62.06           C  
ATOM    665  N9    A A  31      75.083  42.164  -0.541  1.00 56.88           N  
ATOM    666  C8    A A  31      74.083  42.881  -1.152  1.00 56.22           C  
ATOM    667  N7    A A  31      72.877  42.454  -0.869  1.00 55.69           N  
ATOM    668  C5    A A  31      73.093  41.379  -0.019  1.00 52.78           C  
ATOM    669  C6    A A  31      72.210  40.498   0.642  1.00 53.56           C  
ATOM    670  N6    A A  31      70.867  40.568   0.536  1.00 51.52           N  
ATOM    671  N1    A A  31      72.757  39.533   1.428  1.00 53.95           N  
ATOM    672  C2    A A  31      74.103  39.473   1.530  1.00 54.22           C  
ATOM    673  N3    A A  31      75.027  40.244   0.958  1.00 53.01           N  
ATOM    674  C4    A A  31      74.453  41.186   0.189  1.00 54.24           C  
HETATM  675  N1  OMC A  32      74.486  38.301  -2.478  1.00 64.48           N  
HETATM  676  C2  OMC A  32      73.347  37.833  -1.784  1.00 63.86           C  
HETATM  677  N3  OMC A  32      72.174  38.492  -1.928  1.00 62.03           N  
HETATM  678  C4  OMC A  32      72.103  39.570  -2.726  1.00 62.91           C  
HETATM  679  C5  OMC A  32      73.242  40.067  -3.425  1.00 63.75           C  
HETATM  680  C6  OMC A  32      74.401  39.413  -3.269  1.00 64.43           C  
HETATM  681  O2  OMC A  32      73.450  36.817  -1.056  1.00 63.74           O  
HETATM  682  N4  OMC A  32      70.914  40.189  -2.874  1.00 60.56           N  
HETATM  683  C1* OMC A  32      75.756  37.560  -2.358  1.00 68.30           C  
HETATM  684  C2* OMC A  32      75.912  36.485  -3.446  1.00 68.41           C  
HETATM  685  O2* OMC A  32      76.521  35.341  -2.874  1.00 68.27           O  
HETATM  686  CM2 OMC A  32      75.491  34.476  -2.439  1.00 68.09           C  
HETATM  687  C3* OMC A  32      76.801  37.196  -4.464  1.00 69.68           C  
HETATM  688  C4* OMC A  32      77.712  38.038  -3.579  1.00 69.87           C  
HETATM  689  O4* OMC A  32      76.843  38.468  -2.493  1.00 69.46           O  
HETATM  690  O3* OMC A  32      77.529  36.305  -5.304  1.00 71.71           O  
HETATM  691  C5* OMC A  32      78.349  39.239  -4.237  1.00 68.92           C  
HETATM  692  O5* OMC A  32      77.359  40.004  -4.930  1.00 69.40           O  
HETATM  693  P   OMC A  32      77.644  41.515  -5.330  1.00 70.88           P  
HETATM  694  O1P OMC A  32      76.360  42.128  -5.750  1.00 69.87           O  
HETATM  695  O2P OMC A  32      78.813  41.581  -6.239  1.00 71.33           O  
ATOM    696  P     U A  33      76.971  35.962  -6.783  1.00 72.91           P  
ATOM    697  O1P   U A  33      78.044  35.151  -7.420  1.00 73.35           O  
ATOM    698  O2P   U A  33      76.437  37.170  -7.490  1.00 71.97           O  
ATOM    699  O5*   U A  33      75.699  35.053  -6.505  1.00 71.09           O  
ATOM    700  C5*   U A  33      75.822  33.741  -5.951  1.00 71.16           C  
ATOM    701  C4*   U A  33      74.457  33.095  -5.870  1.00 72.48           C  
ATOM    702  O4*   U A  33      73.655  33.729  -4.828  1.00 73.47           O  
ATOM    703  C3*   U A  33      73.601  33.208  -7.129  1.00 72.70           C  
ATOM    704  O3*   U A  33      73.935  32.212  -8.110  1.00 70.79           O  
ATOM    705  C2*   U A  33      72.190  33.041  -6.575  1.00 72.85           C  
ATOM    706  O2*   U A  33      71.917  31.668  -6.340  1.00 73.80           O  
ATOM    707  C1*   U A  33      72.289  33.769  -5.226  1.00 72.04           C  
ATOM    708  N1    U A  33      71.847  35.177  -5.291  1.00 69.22           N  
ATOM    709  C2    U A  33      70.504  35.476  -4.990  1.00 67.66           C  
ATOM    710  O2    U A  33      69.696  34.637  -4.613  1.00 66.73           O  
ATOM    711  N3    U A  33      70.155  36.799  -5.139  1.00 64.32           N  
ATOM    712  C4    U A  33      70.975  37.841  -5.535  1.00 64.94           C  
ATOM    713  O4    U A  33      70.494  38.966  -5.720  1.00 60.96           O  
ATOM    714  C5    U A  33      72.339  37.462  -5.787  1.00 65.88           C  
ATOM    715  C6    U A  33      72.718  36.181  -5.660  1.00 68.66           C  
HETATM  716  P   OMG A  34      73.785  32.556  -9.678  1.00 71.12           P  
HETATM  717  O1P OMG A  34      74.725  33.685 -10.014  1.00 67.66           O  
HETATM  718  O2P OMG A  34      73.921  31.249 -10.387  1.00 68.06           O  
HETATM  719  O5* OMG A  34      72.274  33.065  -9.810  1.00 65.26           O  
HETATM  720  C5* OMG A  34      71.764  33.656 -11.016  1.00 64.26           C  
HETATM  721  C4* OMG A  34      70.295  33.326 -11.163  1.00 64.62           C  
HETATM  722  O4* OMG A  34      70.126  31.890 -11.328  1.00 63.80           O  
HETATM  723  C3* OMG A  34      69.463  33.672  -9.939  1.00 65.81           C  
HETATM  724  O3* OMG A  34      69.094  35.051  -9.965  1.00 67.70           O  
HETATM  725  C2* OMG A  34      68.311  32.644  -9.969  1.00 64.60           C  
HETATM  726  O2* OMG A  34      67.083  32.909 -10.666  1.00 63.66           O  
HETATM  727  CM2 OMG A  34      67.294  33.412 -11.976  1.00 64.30           C  
HETATM  728  C1* OMG A  34      68.999  31.433 -10.593  1.00 61.80           C  
HETATM  729  N9  OMG A  34      69.429  30.396  -9.649  1.00 59.16           N  
HETATM  730  C8  OMG A  34      70.718  29.977  -9.388  1.00 57.76           C  
HETATM  731  N7  OMG A  34      70.774  28.948  -8.579  1.00 56.91           N  
HETATM  732  C5  OMG A  34      69.441  28.691  -8.266  1.00 58.76           C  
HETATM  733  C6  OMG A  34      68.866  27.676  -7.438  1.00 60.12           C  
HETATM  734  O6  OMG A  34      69.451  26.785  -6.793  1.00 60.12           O  
HETATM  735  N1  OMG A  34      67.469  27.777  -7.395  1.00 59.19           N  
HETATM  736  C2  OMG A  34      66.720  28.751  -8.040  1.00 59.35           C  
HETATM  737  N2  OMG A  34      65.376  28.719  -7.843  1.00 56.71           N  
HETATM  738  N3  OMG A  34      67.250  29.698  -8.816  1.00 58.06           N  
HETATM  739  C4  OMG A  34      68.602  29.600  -8.889  1.00 58.63           C  
ATOM    740  P     A A  35      69.287  35.959  -8.637  1.00 69.94           P  
ATOM    741  O1P   A A  35      68.963  37.372  -8.988  1.00 66.80           O  
ATOM    742  O2P   A A  35      70.613  35.629  -8.027  1.00 67.49           O  
ATOM    743  O5*   A A  35      68.138  35.409  -7.674  1.00 69.12           O  
ATOM    744  C5*   A A  35      66.791  35.233  -8.158  1.00 71.95           C  
ATOM    745  C4*   A A  35      65.999  34.348  -7.216  1.00 73.29           C  
ATOM    746  O4*   A A  35      66.330  32.957  -7.426  1.00 71.52           O  
ATOM    747  C3*   A A  35      66.228  34.575  -5.733  1.00 75.03           C  
ATOM    748  O3*   A A  35      65.477  35.689  -5.246  1.00 79.27           O  
ATOM    749  C2*   A A  35      65.806  33.241  -5.122  1.00 73.35           C  
ATOM    750  O2*   A A  35      64.415  33.130  -4.902  1.00 73.82           O  
ATOM    751  C1*   A A  35      66.202  32.251  -6.217  1.00 69.82           C  
ATOM    752  N9    A A  35      67.458  31.568  -5.953  1.00 66.95           N  
ATOM    753  C8    A A  35      68.714  31.908  -6.382  1.00 65.90           C  
ATOM    754  N7    A A  35      69.650  31.095  -5.969  1.00 63.58           N  
ATOM    755  C5    A A  35      68.964  30.158  -5.217  1.00 64.51           C  
ATOM    756  C6    A A  35      69.394  29.030  -4.505  1.00 64.64           C  
ATOM    757  N6    A A  35      70.676  28.646  -4.442  1.00 64.41           N  
ATOM    758  N1    A A  35      68.454  28.302  -3.853  1.00 64.67           N  
ATOM    759  C2    A A  35      67.173  28.694  -3.936  1.00 64.93           C  
ATOM    760  N3    A A  35      66.650  29.737  -4.582  1.00 64.37           N  
ATOM    761  C4    A A  35      67.612  30.435  -5.203  1.00 64.99           C  
ATOM    762  P     A A  36      66.113  36.647  -4.115  1.00 81.67           P  
ATOM    763  O1P   A A  36      65.045  37.623  -3.738  1.00 81.59           O  
ATOM    764  O2P   A A  36      67.432  37.156  -4.600  1.00 82.66           O  
ATOM    765  O5*   A A  36      66.374  35.639  -2.911  1.00 83.98           O  
ATOM    766  C5*   A A  36      65.302  34.813  -2.451  1.00 87.71           C  
ATOM    767  C4*   A A  36      65.801  33.698  -1.565  1.00 90.67           C  
ATOM    768  O4*   A A  36      66.568  32.696  -2.288  1.00 89.82           O  
ATOM    769  C3*   A A  36      66.712  33.999  -0.388  1.00 92.78           C  
ATOM    770  O3*   A A  36      66.031  34.644   0.701  1.00 97.07           O  
ATOM    771  C2*   A A  36      67.162  32.587   0.010  1.00 91.11           C  
ATOM    772  O2*   A A  36      66.278  31.932   0.896  1.00 92.17           O  
ATOM    773  C1*   A A  36      67.143  31.833  -1.328  1.00 88.90           C  
ATOM    774  N9    A A  36      68.488  31.414  -1.719  1.00 85.94           N  
ATOM    775  C8    A A  36      69.437  32.000  -2.530  1.00 84.74           C  
ATOM    776  N7    A A  36      70.588  31.361  -2.537  1.00 83.82           N  
ATOM    777  C5    A A  36      70.369  30.267  -1.700  1.00 84.29           C  
ATOM    778  C6    A A  36      71.195  29.192  -1.278  1.00 84.05           C  
ATOM    779  N6    A A  36      72.462  29.017  -1.673  1.00 85.59           N  
ATOM    780  N1    A A  36      70.658  28.283  -0.430  1.00 83.44           N  
ATOM    781  C2    A A  36      69.385  28.429  -0.052  1.00 82.14           C  
ATOM    782  N3    A A  36      68.512  29.376  -0.383  1.00 82.89           N  
ATOM    783  C4    A A  36      69.074  30.277  -1.213  1.00 84.44           C  
HETATM  784  N1   YG A  37      73.405  29.912   1.754  1.00 91.04           N  
HETATM  785  N2   YG A  37      73.444  28.355   3.405  1.00 90.89           N  
HETATM  786  C2   YG A  37      72.723  29.290   2.751  1.00 90.80           C  
HETATM  787  N3   YG A  37      71.471  29.563   3.082  1.00 90.20           N  
HETATM  788  C3   YG A  37      70.758  28.811   4.134  1.00 89.88           C  
HETATM  789  C4   YG A  37      70.959  30.583   2.344  1.00 90.19           C  
HETATM  790  C5   YG A  37      71.580  31.301   1.340  1.00 89.89           C  
HETATM  791  C6   YG A  37      72.913  30.962   0.974  1.00 90.32           C  
HETATM  792  O6   YG A  37      73.618  31.473   0.097  1.00 89.77           O  
HETATM  793  N7   YG A  37      70.739  32.287   0.844  1.00 89.86           N  
HETATM  794  C8   YG A  37      69.640  32.153   1.532  1.00 89.85           C  
HETATM  795  N9   YG A  37      69.698  31.144   2.464  1.00 90.90           N  
HETATM  796  C10  YG A  37      75.897  27.876   3.472  1.00 91.27           C  
HETATM  797  C11  YG A  37      74.717  28.526   2.928  1.00 92.18           C  
HETATM  798  C12  YG A  37      74.703  29.466   1.794  1.00 92.82           C  
HETATM  799  C13  YG A  37      75.894  29.679   0.821  1.00 95.17           C  
HETATM  800  C14  YG A  37      75.821  28.497  -0.184  1.00 98.02           C  
HETATM  801  C15  YG A  37      76.173  28.762  -1.668  1.00 99.84           C  
HETATM  802  C16  YG A  37      75.991  27.380  -2.269  1.00 99.97           C  
HETATM  803  O17  YG A  37      74.838  26.920  -2.376  1.00100.19           O  
HETATM  804  O18  YG A  37      76.976  26.693  -2.657  1.00100.19           O  
HETATM  805  C19  YG A  37      77.769  27.170  -3.764  1.00 99.79           C  
HETATM  806  N20  YG A  37      75.234  29.792  -2.267  1.00100.19           N  
HETATM  807  C21  YG A  37      75.173  30.151  -3.610  1.00100.19           C  
HETATM  808  O22  YG A  37      74.112  30.084  -4.264  1.00100.19           O  
HETATM  809  O23  YG A  37      76.221  30.547  -4.170  1.00100.19           O  
HETATM  810  C24  YG A  37      76.863  31.711  -3.637  1.00 99.34           C  
HETATM  811  C1*  YG A  37      68.645  30.761   3.410  1.00 93.08           C  
HETATM  812  C2*  YG A  37      68.983  31.220   4.830  1.00 94.21           C  
HETATM  813  O2*  YG A  37      68.545  30.264   5.774  1.00 94.06           O  
HETATM  814  C3*  YG A  37      68.238  32.545   4.908  1.00 94.74           C  
HETATM  815  O3*  YG A  37      68.016  32.978   6.242  1.00 95.58           O  
HETATM  816  C4*  YG A  37      66.976  32.279   4.107  1.00 94.30           C  
HETATM  817  O4*  YG A  37      67.424  31.399   3.039  1.00 94.11           O  
HETATM  818  C5*  YG A  37      66.323  33.509   3.519  1.00 95.27           C  
HETATM  819  O5*  YG A  37      67.311  34.330   2.856  1.00 96.96           O  
HETATM  820  P    YG A  37      66.875  35.460   1.815  1.00 99.01           P  
HETATM  821  O1P  YG A  37      68.123  35.981   1.180  1.00 99.22           O  
HETATM  822  O2P  YG A  37      65.932  36.421   2.473  1.00 99.27           O  
ATOM    823  P     A A  38      69.075  33.984   6.915  1.00 95.95           P  
ATOM    824  O1P   A A  38      68.607  34.325   8.279  1.00 96.86           O  
ATOM    825  O2P   A A  38      69.344  35.071   5.940  1.00 96.06           O  
ATOM    826  O5*   A A  38      70.393  33.112   7.051  1.00 93.01           O  
ATOM    827  C5*   A A  38      70.417  31.929   7.859  1.00 91.33           C  
ATOM    828  C4*   A A  38      71.783  31.310   7.804  1.00 90.73           C  
ATOM    829  O4*   A A  38      72.019  30.786   6.472  1.00 89.61           O  
ATOM    830  C3*   A A  38      72.876  32.339   8.008  1.00 90.95           C  
ATOM    831  O3*   A A  38      73.107  32.623   9.369  1.00 92.09           O  
ATOM    832  C2*   A A  38      74.046  31.768   7.218  1.00 90.07           C  
ATOM    833  O2*   A A  38      74.791  30.775   7.903  1.00 89.63           O  
ATOM    834  C1*   A A  38      73.308  31.178   6.015  1.00 88.20           C  
ATOM    835  N9    A A  38      73.101  32.188   4.971  1.00 84.93           N  
ATOM    836  C8    A A  38      71.911  32.814   4.660  1.00 83.84           C  
ATOM    837  N7    A A  38      72.004  33.675   3.679  1.00 82.59           N  
ATOM    838  C5    A A  38      73.344  33.617   3.312  1.00 82.67           C  
ATOM    839  C6    A A  38      74.081  34.292   2.315  1.00 81.67           C  
ATOM    840  N6    A A  38      73.539  35.181   1.471  1.00 80.44           N  
ATOM    841  N1    A A  38      75.404  34.015   2.214  1.00 80.84           N  
ATOM    842  C2    A A  38      75.938  33.108   3.055  1.00 81.31           C  
ATOM    843  N3    A A  38      75.346  32.402   4.026  1.00 82.09           N  
ATOM    844  C4    A A  38      74.035  32.708   4.104  1.00 83.08           C  
HETATM  845  N1  PSU A  39      74.080  36.066   5.459  1.00 75.82           N  
HETATM  846  C2  PSU A  39      74.415  36.835   4.354  1.00 75.59           C  
HETATM  847  N3  PSU A  39      75.735  36.769   3.984  1.00 76.29           N  
HETATM  848  C4  PSU A  39      76.728  36.038   4.591  1.00 77.28           C  
HETATM  849  C5  PSU A  39      76.307  35.280   5.732  1.00 77.93           C  
HETATM  850  C6  PSU A  39      75.025  35.316   6.112  1.00 76.07           C  
HETATM  851  O2  PSU A  39      73.605  37.525   3.749  1.00 75.80           O  
HETATM  852  O4  PSU A  39      77.875  36.079   4.134  1.00 77.81           O  
HETATM  853  C1* PSU A  39      77.325  34.455   6.488  1.00 79.85           C  
HETATM  854  C2* PSU A  39      78.240  35.315   7.366  1.00 80.79           C  
HETATM  855  O2* PSU A  39      79.550  34.775   7.399  1.00 79.82           O  
HETATM  856  C3* PSU A  39      77.509  35.235   8.700  1.00 81.32           C  
HETATM  857  C4* PSU A  39      77.034  33.800   8.726  1.00 83.30           C  
HETATM  858  O3* PSU A  39      78.312  35.525   9.823  1.00 80.38           O  
HETATM  859  O4* PSU A  39      76.648  33.545   7.349  1.00 81.10           O  
HETATM  860  C5* PSU A  39      75.867  33.557   9.646  1.00 86.38           C  
HETATM  861  O5* PSU A  39      74.796  34.470   9.341  1.00 90.31           O  
HETATM  862  P   PSU A  39      73.308  34.151   9.814  1.00 92.76           P  
HETATM  863  O1P PSU A  39      73.270  34.200  11.303  1.00 93.12           O  
HETATM  864  O2P PSU A  39      72.370  34.998   9.024  1.00 92.11           O  
HETATM  865  P   5MC A  40      78.152  36.942  10.553  1.00 79.24           P  
HETATM  866  O1P 5MC A  40      76.685  37.167  10.693  1.00 80.02           O  
HETATM  867  O2P 5MC A  40      79.021  36.942  11.759  1.00 78.31           O  
HETATM  868  O5* 5MC A  40      78.720  37.994   9.494  1.00 76.52           O  
HETATM  869  C5* 5MC A  40      80.116  38.022   9.141  1.00 72.60           C  
HETATM  870  C4* 5MC A  40      80.351  38.970   7.985  1.00 69.92           C  
HETATM  871  O4* 5MC A  40      79.612  38.525   6.814  1.00 68.54           O  
HETATM  872  C3* 5MC A  40      79.877  40.397   8.200  1.00 68.47           C  
HETATM  873  O3* 5MC A  40      80.825  41.181   8.913  1.00 66.73           O  
HETATM  874  C2* 5MC A  40      79.663  40.906   6.778  1.00 69.16           C  
HETATM  875  O2* 5MC A  40      80.841  41.385   6.165  1.00 68.77           O  
HETATM  876  C1* 5MC A  40      79.168  39.648   6.065  1.00 67.58           C  
HETATM  877  N1  5MC A  40      77.695  39.634   6.019  1.00 65.59           N  
HETATM  878  C2  5MC A  40      77.041  40.408   5.050  1.00 63.72           C  
HETATM  879  O2  5MC A  40      77.718  41.074   4.247  1.00 61.25           O  
HETATM  880  N3  5MC A  40      75.695  40.421   5.017  1.00 63.28           N  
HETATM  881  C4  5MC A  40      74.995  39.705   5.898  1.00 62.84           C  
HETATM  882  N4  5MC A  40      73.662  39.740   5.803  1.00 62.95           N  
HETATM  883  C5  5MC A  40      75.630  38.920   6.903  1.00 63.44           C  
HETATM  884  C6  5MC A  40      76.968  38.900   6.919  1.00 64.81           C  
HETATM  885  CM5 5MC A  40      74.837  38.153   7.923  1.00 63.46           C  
ATOM    886  P     U A  41      80.313  42.447   9.767  1.00 63.29           P  
ATOM    887  O1P   U A  41      81.427  42.965  10.609  1.00 65.66           O  
ATOM    888  O2P   U A  41      79.049  42.016  10.412  1.00 63.65           O  
ATOM    889  O5*   U A  41      79.984  43.552   8.669  1.00 62.24           O  
ATOM    890  C5*   U A  41      81.036  44.143   7.918  1.00 56.79           C  
ATOM    891  C4*   U A  41      80.498  45.168   6.953  1.00 56.45           C  
ATOM    892  O4*   U A  41      79.680  44.526   5.942  1.00 53.37           O  
ATOM    893  C3*   U A  41      79.592  46.251   7.508  1.00 54.69           C  
ATOM    894  O3*   U A  41      80.283  47.292   8.160  1.00 55.80           O  
ATOM    895  C2*   U A  41      78.887  46.737   6.251  1.00 52.98           C  
ATOM    896  O2*   U A  41      79.686  47.539   5.407  1.00 54.28           O  
ATOM    897  C1*   U A  41      78.651  45.424   5.527  1.00 49.85           C  
ATOM    898  N1    U A  41      77.358  44.851   5.928  1.00 45.19           N  
ATOM    899  C2    U A  41      76.197  45.430   5.430  1.00 41.52           C  
ATOM    900  O2    U A  41      76.190  46.418   4.723  1.00 42.41           O  
ATOM    901  N3    U A  41      75.046  44.808   5.803  1.00 41.65           N  
ATOM    902  C4    U A  41      74.928  43.710   6.596  1.00 40.75           C  
ATOM    903  O4    U A  41      73.826  43.203   6.735  1.00 49.24           O  
ATOM    904  C5    U A  41      76.154  43.186   7.099  1.00 43.67           C  
ATOM    905  C6    U A  41      77.302  43.767   6.753  1.00 41.90           C  
ATOM    906  P     G A  42      79.705  47.861   9.545  1.00 53.54           P  
ATOM    907  O1P   G A  42      80.846  48.577  10.148  1.00 56.88           O  
ATOM    908  O2P   G A  42      79.060  46.769  10.292  1.00 47.65           O  
ATOM    909  O5*   G A  42      78.638  48.930   9.059  1.00 51.64           O  
ATOM    910  C5*   G A  42      79.001  49.877   8.054  1.00 51.41           C  
ATOM    911  C4*   G A  42      77.769  50.433   7.372  1.00 54.00           C  
ATOM    912  O4*   G A  42      77.032  49.365   6.704  1.00 51.14           O  
ATOM    913  C3*   G A  42      76.726  51.067   8.272  1.00 53.82           C  
ATOM    914  O3*   G A  42      77.060  52.389   8.660  1.00 57.90           O  
ATOM    915  C2*   G A  42      75.480  51.029   7.397  1.00 53.12           C  
ATOM    916  O2*   G A  42      75.444  52.122   6.489  1.00 55.26           O  
ATOM    917  C1*   G A  42      75.657  49.688   6.671  1.00 49.13           C  
ATOM    918  N9    G A  42      74.906  48.601   7.302  1.00 44.46           N  
ATOM    919  C8    G A  42      75.381  47.535   8.030  1.00 41.11           C  
ATOM    920  N7    G A  42      74.431  46.730   8.433  1.00 41.58           N  
ATOM    921  C5    G A  42      73.265  47.311   7.951  1.00 42.30           C  
ATOM    922  C6    G A  42      71.901  46.887   8.049  1.00 39.56           C  
ATOM    923  O6    G A  42      71.457  45.892   8.572  1.00 41.46           O  
ATOM    924  N1    G A  42      71.044  47.776   7.421  1.00 40.73           N  
ATOM    925  C2    G A  42      71.439  48.906   6.760  1.00 39.14           C  
ATOM    926  N2    G A  42      70.481  49.648   6.204  1.00 43.75           N  
ATOM    927  N3    G A  42      72.696  49.297   6.643  1.00 42.69           N  
ATOM    928  C4    G A  42      73.544  48.463   7.263  1.00 41.15           C  
ATOM    929  P     G A  43      76.348  53.024   9.953  1.00 59.25           P  
ATOM    930  O1P   G A  43      76.997  54.350  10.246  1.00 56.83           O  
ATOM    931  O2P   G A  43      76.330  51.946  10.985  1.00 53.33           O  
ATOM    932  O5*   G A  43      74.844  53.245   9.472  1.00 55.31           O  
ATOM    933  C5*   G A  43      74.543  54.211   8.454  1.00 54.14           C  
ATOM    934  C4*   G A  43      73.059  54.264   8.204  1.00 51.95           C  
ATOM    935  O4*   G A  43      72.595  52.976   7.736  1.00 49.20           O  
ATOM    936  C3*   G A  43      72.190  54.541   9.414  1.00 53.09           C  
ATOM    937  O3*   G A  43      72.143  55.927   9.696  1.00 58.30           O  
ATOM    938  C2*   G A  43      70.852  53.982   8.963  1.00 50.69           C  
ATOM    939  O2*   G A  43      70.277  54.842   8.008  1.00 55.19           O  
ATOM    940  C1*   G A  43      71.295  52.730   8.219  1.00 46.94           C  
ATOM    941  N9    G A  43      71.350  51.546   9.071  1.00 44.87           N  
ATOM    942  C8    G A  43      72.451  50.993   9.658  1.00 41.57           C  
ATOM    943  N7    G A  43      72.182  49.898  10.321  1.00 41.08           N  
ATOM    944  C5    G A  43      70.826  49.729  10.161  1.00 40.57           C  
ATOM    945  C6    G A  43      69.965  48.713  10.615  1.00 40.78           C  
ATOM    946  O6    G A  43      70.230  47.688  11.221  1.00 38.22           O  
ATOM    947  N1    G A  43      68.652  48.964  10.259  1.00 41.30           N  
ATOM    948  C2    G A  43      68.217  50.037   9.555  1.00 42.51           C  
ATOM    949  N2    G A  43      66.869  50.119   9.354  1.00 42.42           N  
ATOM    950  N3    G A  43      69.008  50.969   9.090  1.00 44.65           N  
ATOM    951  C4    G A  43      70.291  50.756   9.425  1.00 44.20           C  
ATOM    952  P     A A  44      71.830  56.434  11.191  1.00 59.72           P  
ATOM    953  O1P   A A  44      72.143  57.887  11.179  1.00 60.18           O  
ATOM    954  O2P   A A  44      72.477  55.542  12.203  1.00 59.03           O  
ATOM    955  O5*   A A  44      70.251  56.220  11.341  1.00 61.19           O  
ATOM    956  C5*   A A  44      69.317  56.879  10.465  1.00 59.15           C  
ATOM    957  C4*   A A  44      67.971  56.191  10.523  1.00 58.59           C  
ATOM    958  O4*   A A  44      68.135  54.794  10.173  1.00 58.70           O  
ATOM    959  C3*   A A  44      67.277  56.140  11.871  1.00 58.67           C  
ATOM    960  O3*   A A  44      66.523  57.316  12.123  1.00 60.86           O  
ATOM    961  C2*   A A  44      66.365  54.931  11.744  1.00 56.86           C  
ATOM    962  O2*   A A  44      65.162  55.206  11.061  1.00 56.94           O  
ATOM    963  C1*   A A  44      67.208  53.996  10.887  1.00 54.14           C  
ATOM    964  N9    A A  44      67.970  53.042  11.685  1.00 47.88           N  
ATOM    965  C8    A A  44      69.290  53.122  12.063  1.00 44.83           C  
ATOM    966  N7    A A  44      69.697  52.091  12.767  1.00 43.37           N  
ATOM    967  C5    A A  44      68.570  51.276  12.856  1.00 41.45           C  
ATOM    968  C6    A A  44      68.345  50.029  13.433  1.00 42.09           C  
ATOM    969  N6    A A  44      69.278  49.331  14.092  1.00 45.08           N  
ATOM    970  N1    A A  44      67.101  49.482  13.313  1.00 42.78           N  
ATOM    971  C2    A A  44      66.176  50.170  12.653  1.00 41.91           C  
ATOM    972  N3    A A  44      66.268  51.351  12.075  1.00 40.34           N  
ATOM    973  C4    A A  44      67.500  51.858  12.203  1.00 45.75           C  
ATOM    974  P     G A  45      66.470  57.914  13.620  1.00 58.58           P  
ATOM    975  O1P   G A  45      65.444  58.996  13.659  1.00 60.46           O  
ATOM    976  O2P   G A  45      67.862  58.228  14.011  1.00 57.19           O  
ATOM    977  O5*   G A  45      65.905  56.683  14.444  1.00 57.67           O  
ATOM    978  C5*   G A  45      64.533  56.287  14.311  1.00 52.30           C  
ATOM    979  C4*   G A  45      64.255  55.158  15.248  1.00 50.49           C  
ATOM    980  O4*   G A  45      65.017  54.002  14.829  1.00 50.44           O  
ATOM    981  C3*   G A  45      64.685  55.408  16.681  1.00 49.36           C  
ATOM    982  O3*   G A  45      63.616  56.059  17.359  1.00 51.42           O  
ATOM    983  C2*   G A  45      64.894  53.991  17.190  1.00 48.83           C  
ATOM    984  O2*   G A  45      63.645  53.382  17.453  1.00 45.60           O  
ATOM    985  C1*   G A  45      65.456  53.277  15.954  1.00 45.24           C  
ATOM    986  N9    G A  45      66.920  53.202  15.909  1.00 41.15           N  
ATOM    987  C8    G A  45      67.772  54.149  15.386  1.00 40.45           C  
ATOM    988  N7    G A  45      69.031  53.833  15.507  1.00 39.19           N  
ATOM    989  C5    G A  45      69.014  52.601  16.134  1.00 35.16           C  
ATOM    990  C6    G A  45      70.070  51.816  16.551  1.00 38.09           C  
ATOM    991  O6    G A  45      71.292  52.028  16.410  1.00 39.62           O  
ATOM    992  N1    G A  45      69.636  50.685  17.196  1.00 37.05           N  
ATOM    993  C2    G A  45      68.338  50.332  17.390  1.00 38.47           C  
ATOM    994  N2    G A  45      68.164  49.150  17.996  1.00 39.48           N  
ATOM    995  N3    G A  45      67.301  51.072  17.004  1.00 37.37           N  
ATOM    996  C4    G A  45      67.719  52.193  16.390  1.00 36.46           C  
HETATM  997  P   7MG A  46      63.905  57.326  18.310  1.00 53.01           P  
HETATM  998  O1P 7MG A  46      64.951  58.232  17.766  1.00 52.60           O  
HETATM  999  O2P 7MG A  46      62.558  57.883  18.619  1.00 54.84           O  
HETATM 1000  O5* 7MG A  46      64.457  56.673  19.662  1.00 50.92           O  
HETATM 1001  C5* 7MG A  46      63.673  55.697  20.380  1.00 47.09           C  
HETATM 1002  C4* 7MG A  46      63.911  55.833  21.859  1.00 45.87           C  
HETATM 1003  O4* 7MG A  46      65.276  55.388  22.160  1.00 44.37           O  
HETATM 1004  C3* 7MG A  46      63.810  57.269  22.378  1.00 44.64           C  
HETATM 1005  O3* 7MG A  46      63.177  57.259  23.652  1.00 47.38           O  
HETATM 1006  C2* 7MG A  46      65.280  57.700  22.484  1.00 44.78           C  
HETATM 1007  O2* 7MG A  46      65.577  58.744  23.387  1.00 43.73           O  
HETATM 1008  C1* 7MG A  46      65.949  56.384  22.900  1.00 39.69           C  
HETATM 1009  N9  7MG A  46      67.386  56.266  22.628  1.00 35.94           N  
HETATM 1010  C8  7MG A  46      68.201  57.146  21.945  1.00 35.57           C  
HETATM 1011  N7  7MG A  46      69.429  56.680  21.823  1.00 35.71           N  
HETATM 1012  C5  7MG A  46      69.423  55.456  22.475  1.00 30.85           C  
HETATM 1013  C6  7MG A  46      70.472  54.529  22.704  1.00 32.44           C  
HETATM 1014  O6  7MG A  46      71.653  54.665  22.404  1.00 34.18           O  
HETATM 1015  N1  7MG A  46      70.035  53.405  23.371  1.00 27.06           N  
HETATM 1016  C2  7MG A  46      68.726  53.216  23.821  1.00 28.72           C  
HETATM 1017  N2  7MG A  46      68.413  52.076  24.386  1.00 30.51           N  
HETATM 1018  N3  7MG A  46      67.782  54.113  23.692  1.00 29.24           N  
HETATM 1019  C4  7MG A  46      68.183  55.188  22.989  1.00 30.11           C  
HETATM 1020  CM7 7MG A  46      70.529  57.362  21.130  1.00 35.54           C  
ATOM   1021  P     U A  47      61.607  57.698  23.777  1.00 47.19           P  
ATOM   1022  O1P   U A  47      60.804  56.464  23.970  1.00 48.54           O  
ATOM   1023  O2P   U A  47      61.264  58.639  22.661  1.00 50.14           O  
ATOM   1024  O5*   U A  47      61.587  58.598  25.089  1.00 45.96           O  
ATOM   1025  C5*   U A  47      62.341  59.800  25.124  1.00 48.85           C  
ATOM   1026  C4*   U A  47      62.482  60.255  26.545  1.00 50.43           C  
ATOM   1027  O4*   U A  47      61.241  60.814  27.040  1.00 50.06           O  
ATOM   1028  C3*   U A  47      62.833  59.116  27.493  1.00 49.46           C  
ATOM   1029  O3*   U A  47      63.658  59.705  28.465  1.00 50.38           O  
ATOM   1030  C2*   U A  47      61.488  58.710  28.088  1.00 48.84           C  
ATOM   1031  O2*   U A  47      61.620  58.065  29.342  1.00 41.07           O  
ATOM   1032  C1*   U A  47      60.798  60.075  28.167  1.00 50.55           C  
ATOM   1033  N1    U A  47      59.332  60.087  28.149  1.00 54.98           N  
ATOM   1034  C2    U A  47      58.686  60.797  29.157  1.00 56.81           C  
ATOM   1035  O2    U A  47      59.281  61.400  30.035  1.00 56.24           O  
ATOM   1036  N3    U A  47      57.319  60.784  29.079  1.00 58.14           N  
ATOM   1037  C4    U A  47      56.552  60.165  28.123  1.00 58.07           C  
ATOM   1038  O4    U A  47      55.332  60.325  28.140  1.00 59.44           O  
ATOM   1039  C5    U A  47      57.289  59.452  27.121  1.00 57.49           C  
ATOM   1040  C6    U A  47      58.613  59.436  27.164  1.00 55.97           C  
ATOM   1041  P     C A  48      65.235  59.648  28.260  1.00 51.91           P  
ATOM   1042  O1P   C A  48      65.823  60.627  29.208  1.00 48.27           O  
ATOM   1043  O2P   C A  48      65.455  59.769  26.780  1.00 53.23           O  
ATOM   1044  O5*   C A  48      65.585  58.155  28.696  1.00 48.16           O  
ATOM   1045  C5*   C A  48      65.553  57.821  30.068  1.00 37.25           C  
ATOM   1046  C4*   C A  48      66.046  56.424  30.251  1.00 33.96           C  
ATOM   1047  O4*   C A  48      67.242  56.203  29.430  1.00 33.29           O  
ATOM   1048  C3*   C A  48      66.457  56.168  31.682  1.00 32.27           C  
ATOM   1049  O3*   C A  48      66.149  54.821  32.011  1.00 29.89           O  
ATOM   1050  C2*   C A  48      67.972  56.366  31.634  1.00 27.98           C  
ATOM   1051  O2*   C A  48      68.646  55.696  32.665  1.00 27.27           O  
ATOM   1052  C1*   C A  48      68.291  55.799  30.245  1.00 27.98           C  
ATOM   1053  N1    C A  48      69.554  56.279  29.660  1.00 23.61           N  
ATOM   1054  C2    C A  48      70.690  55.421  29.693  1.00 24.90           C  
ATOM   1055  O2    C A  48      70.567  54.326  30.181  1.00 24.38           O  
ATOM   1056  N3    C A  48      71.884  55.875  29.215  1.00 23.05           N  
ATOM   1057  C4    C A  48      71.964  57.127  28.722  1.00 24.65           C  
ATOM   1058  N4    C A  48      73.163  57.592  28.288  1.00 21.79           N  
ATOM   1059  C5    C A  48      70.843  57.978  28.653  1.00 27.10           C  
ATOM   1060  C6    C A  48      69.665  57.518  29.131  1.00 23.05           C  
HETATM 1061  P   5MC A  49      65.638  54.464  33.461  1.00 31.61           P  
HETATM 1062  O1P 5MC A  49      65.643  52.989  33.540  1.00 33.43           O  
HETATM 1063  O2P 5MC A  49      66.252  55.251  34.624  1.00 26.54           O  
HETATM 1064  O5* 5MC A  49      64.126  54.972  33.441  1.00 31.88           O  
HETATM 1065  C5* 5MC A  49      63.204  54.605  32.392  1.00 32.68           C  
HETATM 1066  C4* 5MC A  49      61.796  55.006  32.810  1.00 30.80           C  
HETATM 1067  O4* 5MC A  49      61.292  54.110  33.848  1.00 32.40           O  
HETATM 1068  C3* 5MC A  49      61.745  56.381  33.468  1.00 31.44           C  
HETATM 1069  O3* 5MC A  49      61.666  57.398  32.483  1.00 31.25           O  
HETATM 1070  C2* 5MC A  49      60.485  56.287  34.312  1.00 34.62           C  
HETATM 1071  O2* 5MC A  49      59.414  56.469  33.417  1.00 35.41           O  
HETATM 1072  C1* 5MC A  49      60.568  54.849  34.825  1.00 33.10           C  
HETATM 1073  N1  5MC A  49      61.331  54.787  36.102  1.00 32.34           N  
HETATM 1074  C2  5MC A  49      60.737  55.340  37.271  1.00 32.42           C  
HETATM 1075  O2  5MC A  49      59.568  55.832  37.197  1.00 28.24           O  
HETATM 1076  N3  5MC A  49      61.420  55.348  38.428  1.00 29.50           N  
HETATM 1077  C4  5MC A  49      62.626  54.785  38.503  1.00 32.97           C  
HETATM 1078  N4  5MC A  49      63.188  54.744  39.702  1.00 28.75           N  
HETATM 1079  C5  5MC A  49      63.270  54.215  37.358  1.00 34.81           C  
HETATM 1080  C6  5MC A  49      62.583  54.243  36.161  1.00 33.69           C  
HETATM 1081  CM5 5MC A  49      64.655  53.581  37.465  1.00 34.52           C  
ATOM   1082  P     U A  50      62.421  58.828  32.744  1.00 33.43           P  
ATOM   1083  O1P   U A  50      62.480  59.502  31.468  1.00 33.78           O  
ATOM   1084  O2P   U A  50      63.651  58.683  33.637  1.00 30.15           O  
ATOM   1085  O5*   U A  50      61.374  59.605  33.659  1.00 33.34           O  
ATOM   1086  C5*   U A  50      60.011  59.728  33.268  1.00 34.47           C  
ATOM   1087  C4*   U A  50      59.218  60.444  34.331  1.00 31.72           C  
ATOM   1088  O4*   U A  50      58.928  59.552  35.438  1.00 29.81           O  
ATOM   1089  C3*   U A  50      59.871  61.655  34.983  1.00 32.89           C  
ATOM   1090  O3*   U A  50      59.731  62.817  34.131  1.00 33.31           O  
ATOM   1091  C2*   U A  50      59.066  61.725  36.279  1.00 33.15           C  
ATOM   1092  O2*   U A  50      57.727  62.115  35.997  1.00 35.03           O  
ATOM   1093  C1*   U A  50      58.966  60.259  36.646  1.00 31.45           C  
ATOM   1094  N1    U A  50      60.113  59.812  37.453  1.00 31.52           N  
ATOM   1095  C2    U A  50      60.106  60.183  38.793  1.00 30.82           C  
ATOM   1096  O2    U A  50      59.256  60.941  39.253  1.00 30.94           O  
ATOM   1097  N3    U A  50      61.128  59.648  39.560  1.00 28.21           N  
ATOM   1098  C4    U A  50      62.143  58.831  39.099  1.00 27.57           C  
ATOM   1099  O4    U A  50      62.938  58.308  39.924  1.00 33.33           O  
ATOM   1100  C5    U A  50      62.120  58.581  37.684  1.00 29.31           C  
ATOM   1101  C6    U A  50      61.140  59.068  36.927  1.00 31.19           C  
ATOM   1102  P     G A  51      60.854  63.992  34.167  1.00 39.01           P  
ATOM   1103  O1P   G A  51      60.476  64.987  33.121  1.00 35.87           O  
ATOM   1104  O2P   G A  51      62.217  63.402  34.134  1.00 37.91           O  
ATOM   1105  O5*   G A  51      60.648  64.595  35.600  1.00 36.58           O  
ATOM   1106  C5*   G A  51      59.391  65.231  35.906  1.00 40.92           C  
ATOM   1107  C4*   G A  51      59.365  65.668  37.334  1.00 38.83           C  
ATOM   1108  O4*   G A  51      59.409  64.501  38.203  1.00 38.03           O  
ATOM   1109  C3*   G A  51      60.529  66.527  37.797  1.00 40.34           C  
ATOM   1110  O3*   G A  51      60.383  67.887  37.400  1.00 41.89           O  
ATOM   1111  C2*   G A  51      60.488  66.292  39.298  1.00 38.39           C  
ATOM   1112  O2*   G A  51      59.369  66.938  39.873  1.00 35.78           O  
ATOM   1113  C1*   G A  51      60.143  64.799  39.378  1.00 37.87           C  
ATOM   1114  N9    G A  51      61.303  63.902  39.460  1.00 32.83           N  
ATOM   1115  C8    G A  51      61.915  63.238  38.424  1.00 33.29           C  
ATOM   1116  N7    G A  51      62.917  62.483  38.830  1.00 36.76           N  
ATOM   1117  C5    G A  51      62.986  62.690  40.208  1.00 35.57           C  
ATOM   1118  C6    G A  51      63.913  62.188  41.197  1.00 34.80           C  
ATOM   1119  O6    G A  51      64.846  61.404  41.049  1.00 35.08           O  
ATOM   1120  N1    G A  51      63.653  62.706  42.466  1.00 36.79           N  
ATOM   1121  C2    G A  51      62.620  63.573  42.759  1.00 37.65           C  
ATOM   1122  N2    G A  51      62.530  63.985  44.032  1.00 35.39           N  
ATOM   1123  N3    G A  51      61.746  64.007  41.868  1.00 36.74           N  
ATOM   1124  C4    G A  51      61.996  63.548  40.619  1.00 35.55           C  
ATOM   1125  P     U A  52      61.706  68.803  37.121  1.00 42.84           P  
ATOM   1126  O1P   U A  52      61.321  70.142  36.563  1.00 41.91           O  
ATOM   1127  O2P   U A  52      62.726  68.025  36.369  1.00 41.76           O  
ATOM   1128  O5*   U A  52      62.224  68.978  38.598  1.00 36.07           O  
ATOM   1129  C5*   U A  52      61.477  69.672  39.558  1.00 34.35           C  
ATOM   1130  C4*   U A  52      62.118  69.559  40.905  1.00 33.37           C  
ATOM   1131  O4*   U A  52      62.012  68.190  41.401  1.00 32.65           O  
ATOM   1132  C3*   U A  52      63.616  69.845  41.000  1.00 33.39           C  
ATOM   1133  O3*   U A  52      63.932  71.237  41.043  1.00 39.74           O  
ATOM   1134  C2*   U A  52      63.939  69.176  42.335  1.00 30.75           C  
ATOM   1135  O2*   U A  52      63.396  69.967  43.391  1.00 32.96           O  
ATOM   1136  C1*   U A  52      63.119  67.888  42.237  1.00 31.31           C  
ATOM   1137  N1    U A  52      63.946  66.845  41.579  1.00 30.44           N  
ATOM   1138  C2    U A  52      64.863  66.182  42.369  1.00 30.06           C  
ATOM   1139  O2    U A  52      65.009  66.454  43.550  1.00 29.26           O  
ATOM   1140  N3    U A  52      65.623  65.231  41.729  1.00 29.70           N  
ATOM   1141  C4    U A  52      65.587  64.922  40.408  1.00 34.08           C  
ATOM   1142  O4    U A  52      66.347  64.036  39.967  1.00 34.26           O  
ATOM   1143  C5    U A  52      64.622  65.680  39.641  1.00 33.07           C  
ATOM   1144  C6    U A  52      63.831  66.587  40.266  1.00 26.53           C  
ATOM   1145  P     G A  53      65.414  71.746  40.572  1.00 37.07           P  
ATOM   1146  O1P   G A  53      65.308  73.211  40.502  1.00 37.71           O  
ATOM   1147  O2P   G A  53      65.813  70.998  39.374  1.00 30.94           O  
ATOM   1148  O5*   G A  53      66.393  71.211  41.694  1.00 32.28           O  
ATOM   1149  C5*   G A  53      66.229  71.538  43.081  1.00 34.99           C  
ATOM   1150  C4*   G A  53      67.174  70.724  43.912  1.00 34.42           C  
ATOM   1151  O4*   G A  53      66.864  69.323  43.765  1.00 35.23           O  
ATOM   1152  C3*   G A  53      68.665  70.782  43.620  1.00 34.67           C  
ATOM   1153  O3*   G A  53      69.247  71.895  44.262  1.00 38.95           O  
ATOM   1154  C2*   G A  53      69.164  69.500  44.283  1.00 34.61           C  
ATOM   1155  O2*   G A  53      69.182  69.612  45.688  1.00 34.19           O  
ATOM   1156  C1*   G A  53      68.019  68.537  43.991  1.00 32.96           C  
ATOM   1157  N9    G A  53      68.320  67.761  42.786  1.00 28.24           N  
ATOM   1158  C8    G A  53      67.793  67.893  41.556  1.00 31.68           C  
ATOM   1159  N7    G A  53      68.284  67.042  40.701  1.00 26.28           N  
ATOM   1160  C5    G A  53      69.196  66.319  41.426  1.00 28.92           C  
ATOM   1161  C6    G A  53      70.091  65.297  41.013  1.00 25.46           C  
ATOM   1162  O6    G A  53      70.220  64.815  39.895  1.00 28.07           O  
ATOM   1163  N1    G A  53      70.897  64.869  42.059  1.00 29.47           N  
ATOM   1164  C2    G A  53      70.866  65.396  43.328  1.00 26.11           C  
ATOM   1165  N2    G A  53      71.750  64.920  44.203  1.00 30.39           N  
ATOM   1166  N3    G A  53      70.033  66.332  43.705  1.00 30.28           N  
ATOM   1167  C4    G A  53      69.234  66.743  42.708  1.00 26.13           C  
HETATM 1168  N1  5MU A  54      73.251  67.803  42.694  1.00 35.32           N  
HETATM 1169  C2  5MU A  54      73.839  66.792  41.983  1.00 34.56           C  
HETATM 1170  N3  5MU A  54      73.442  66.712  40.672  1.00 35.04           N  
HETATM 1171  C4  5MU A  54      72.524  67.496  40.026  1.00 33.29           C  
HETATM 1172  C5  5MU A  54      71.922  68.529  40.835  1.00 31.37           C  
HETATM 1173  C5M 5MU A  54      70.918  69.466  40.203  1.00 28.02           C  
HETATM 1174  C6  5MU A  54      72.301  68.618  42.113  1.00 33.12           C  
HETATM 1175  O2  5MU A  54      74.649  66.016  42.473  1.00 42.24           O  
HETATM 1176  O4  5MU A  54      72.257  67.283  38.828  1.00 36.52           O  
HETATM 1177  C1* 5MU A  54      73.711  68.026  44.082  1.00 37.46           C  
HETATM 1178  C2* 5MU A  54      75.024  68.864  44.109  1.00 40.51           C  
HETATM 1179  O2* 5MU A  54      75.836  68.430  45.180  1.00 39.54           O  
HETATM 1180  C3* 5MU A  54      74.485  70.281  44.321  1.00 42.37           C  
HETATM 1181  C4* 5MU A  54      73.295  70.030  45.246  1.00 39.47           C  
HETATM 1182  O3* 5MU A  54      75.448  71.209  44.879  1.00 44.31           O  
HETATM 1183  O4* 5MU A  54      72.728  68.779  44.751  1.00 36.92           O  
HETATM 1184  C5* 5MU A  54      72.225  71.084  45.312  1.00 38.24           C  
HETATM 1185  O5* 5MU A  54      71.693  71.314  44.026  1.00 36.39           O  
HETATM 1186  P   5MU A  54      70.668  72.485  43.743  1.00 37.08           P  
HETATM 1187  O1P 5MU A  54      70.657  72.754  42.270  1.00 33.03           O  
HETATM 1188  O2P 5MU A  54      70.866  73.636  44.695  1.00 39.97           O  
HETATM 1189  N1  PSU A  55      74.158  70.927  39.519  1.00 35.82           N  
HETATM 1190  C2  PSU A  55      73.717  70.455  38.323  1.00 38.30           C  
HETATM 1191  N3  PSU A  55      74.479  69.441  37.783  1.00 35.34           N  
HETATM 1192  C4  PSU A  55      75.687  68.934  38.291  1.00 36.07           C  
HETATM 1193  C5  PSU A  55      76.107  69.537  39.499  1.00 33.56           C  
HETATM 1194  C6  PSU A  55      75.337  70.458  40.076  1.00 35.52           C  
HETATM 1195  O2  PSU A  55      72.728  70.924  37.738  1.00 37.13           O  
HETATM 1196  O4  PSU A  55      76.304  68.043  37.675  1.00 32.17           O  
HETATM 1197  C1* PSU A  55      77.461  69.118  40.100  1.00 34.24           C  
HETATM 1198  C2* PSU A  55      78.634  70.014  39.665  1.00 37.02           C  
HETATM 1199  O2* PSU A  55      79.793  69.181  39.668  1.00 37.99           O  
HETATM 1200  C3* PSU A  55      78.650  71.033  40.796  1.00 38.47           C  
HETATM 1201  C4* PSU A  55      78.398  70.137  41.999  1.00 35.54           C  
HETATM 1202  O3* PSU A  55      79.864  71.807  40.930  1.00 37.81           O  
HETATM 1203  O4* PSU A  55      77.424  69.161  41.505  1.00 34.94           O  
HETATM 1204  C5* PSU A  55      77.870  70.849  43.223  1.00 36.16           C  
HETATM 1205  O5* PSU A  55      76.820  71.796  42.859  1.00 39.91           O  
HETATM 1206  P   PSU A  55      75.925  72.471  43.991  1.00 43.83           P  
HETATM 1207  O1P PSU A  55      74.766  73.172  43.306  1.00 46.32           O  
HETATM 1208  O2P PSU A  55      76.844  73.293  44.834  1.00 45.38           O  
ATOM   1209  P     C A  56      79.893  73.323  40.427  1.00 37.81           P  
ATOM   1210  O1P   C A  56      81.372  73.714  40.748  1.00 37.74           O  
ATOM   1211  O2P   C A  56      78.768  74.121  41.003  1.00 35.29           O  
ATOM   1212  O5*   C A  56      79.825  73.304  38.871  1.00 33.59           O  
ATOM   1213  C5*   C A  56      79.545  74.494  38.135  1.00 33.82           C  
ATOM   1214  C4*   C A  56      79.735  74.239  36.678  1.00 32.83           C  
ATOM   1215  O4*   C A  56      81.120  73.946  36.360  1.00 33.52           O  
ATOM   1216  C3*   C A  56      78.954  73.036  36.142  1.00 30.99           C  
ATOM   1217  O3*   C A  56      77.596  73.405  35.829  1.00 31.77           O  
ATOM   1218  C2*   C A  56      79.739  72.659  34.901  1.00 31.94           C  
ATOM   1219  O2*   C A  56      79.480  73.400  33.725  1.00 29.83           O  
ATOM   1220  C1*   C A  56      81.174  72.955  35.333  1.00 33.49           C  
ATOM   1221  N1    C A  56      81.800  71.739  35.857  1.00 29.51           N  
ATOM   1222  C2    C A  56      82.298  70.818  34.956  1.00 28.07           C  
ATOM   1223  O2    C A  56      82.160  71.013  33.755  1.00 29.82           O  
ATOM   1224  N3    C A  56      82.907  69.716  35.412  1.00 27.72           N  
ATOM   1225  C4    C A  56      83.035  69.503  36.720  1.00 32.68           C  
ATOM   1226  N4    C A  56      83.675  68.399  37.102  1.00 29.93           N  
ATOM   1227  C5    C A  56      82.515  70.397  37.674  1.00 30.58           C  
ATOM   1228  C6    C A  56      81.894  71.503  37.204  1.00 29.09           C  
ATOM   1229  P     G A  57      76.379  72.458  36.349  1.00 32.12           P  
ATOM   1230  O1P   G A  57      75.189  73.231  35.952  1.00 34.27           O  
ATOM   1231  O2P   G A  57      76.558  72.064  37.745  1.00 22.22           O  
ATOM   1232  O5*   G A  57      76.529  71.173  35.451  1.00 28.84           O  
ATOM   1233  C5*   G A  57      76.321  71.233  34.052  1.00 30.38           C  
ATOM   1234  C4*   G A  57      76.829  69.964  33.385  1.00 32.75           C  
ATOM   1235  O4*   G A  57      78.284  69.786  33.566  1.00 34.57           O  
ATOM   1236  C3*   G A  57      76.221  68.672  33.883  1.00 29.11           C  
ATOM   1237  O3*   G A  57      74.937  68.457  33.239  1.00 27.73           O  
ATOM   1238  C2*   G A  57      77.251  67.650  33.423  1.00 28.35           C  
ATOM   1239  O2*   G A  57      77.136  67.492  32.013  1.00 34.71           O  
ATOM   1240  C1*   G A  57      78.571  68.404  33.655  1.00 32.37           C  
ATOM   1241  N9    G A  57      79.153  68.106  34.968  1.00 26.80           N  
ATOM   1242  C8    G A  57      79.113  68.840  36.137  1.00 26.67           C  
ATOM   1243  N7    G A  57      79.792  68.266  37.110  1.00 27.09           N  
ATOM   1244  C5    G A  57      80.280  67.087  36.529  1.00 28.66           C  
ATOM   1245  C6    G A  57      81.103  66.037  37.080  1.00 28.50           C  
ATOM   1246  O6    G A  57      81.588  65.974  38.203  1.00 28.04           O  
ATOM   1247  N1    G A  57      81.358  65.033  36.153  1.00 27.23           N  
ATOM   1248  C2    G A  57      80.910  65.036  34.858  1.00 27.75           C  
ATOM   1249  N2    G A  57      81.215  63.982  34.112  1.00 25.22           N  
ATOM   1250  N3    G A  57      80.185  66.017  34.328  1.00 27.62           N  
ATOM   1251  C4    G A  57      79.901  66.989  35.227  1.00 25.05           C  
HETATM 1252  P   1MA A  58      73.770  67.765  34.057  1.00 30.65           P  
HETATM 1253  O1P 1MA A  58      73.621  68.229  35.450  1.00 29.49           O  
HETATM 1254  O2P 1MA A  58      72.638  67.886  33.105  1.00 32.84           O  
HETATM 1255  O5* 1MA A  58      74.315  66.273  34.254  1.00 28.81           O  
HETATM 1256  C5* 1MA A  58      74.592  65.439  33.080  1.00 29.42           C  
HETATM 1257  C4* 1MA A  58      74.279  63.972  33.383  1.00 33.42           C  
HETATM 1258  O4* 1MA A  58      74.880  63.685  34.667  1.00 32.36           O  
HETATM 1259  C3* 1MA A  58      72.789  63.573  33.509  1.00 35.13           C  
HETATM 1260  O3* 1MA A  58      72.625  62.168  33.250  1.00 36.80           O  
HETATM 1261  C2* 1MA A  58      72.560  63.667  35.012  1.00 34.80           C  
HETATM 1262  O2* 1MA A  58      71.525  62.828  35.506  1.00 36.27           O  
HETATM 1263  C1* 1MA A  58      73.908  63.150  35.551  1.00 33.62           C  
HETATM 1264  N9  1MA A  58      74.284  63.494  36.930  1.00 30.36           N  
HETATM 1265  C8  1MA A  58      73.887  64.574  37.688  1.00 34.55           C  
HETATM 1266  N7  1MA A  58      74.415  64.610  38.899  1.00 33.32           N  
HETATM 1267  C5  1MA A  58      75.204  63.469  38.953  1.00 33.37           C  
HETATM 1268  C6  1MA A  58      76.031  62.941  39.948  1.00 33.58           C  
HETATM 1269  N6  1MA A  58      76.184  63.488  41.134  1.00 41.19           N  
HETATM 1270  N1  1MA A  58      76.708  61.803  39.669  1.00 34.48           N  
HETATM 1271  CM1 1MA A  58      77.649  61.222  40.626  1.00 31.43           C  
HETATM 1272  C2  1MA A  58      76.527  61.216  38.479  1.00 28.43           C  
HETATM 1273  N3  1MA A  58      75.793  61.624  37.453  1.00 31.67           N  
HETATM 1274  C4  1MA A  58      75.142  62.771  37.747  1.00 33.02           C  
ATOM   1275  P     U A  59      72.617  61.530  31.733  1.00 41.00           P  
ATOM   1276  O1P   U A  59      73.971  61.410  31.109  1.00 33.02           O  
ATOM   1277  O2P   U A  59      71.557  62.222  31.005  1.00 41.76           O  
ATOM   1278  O5*   U A  59      72.130  60.048  31.994  1.00 39.14           O  
ATOM   1279  C5*   U A  59      70.719  59.794  32.346  1.00 35.32           C  
ATOM   1280  C4*   U A  59      70.618  58.472  33.026  1.00 29.87           C  
ATOM   1281  O4*   U A  59      71.242  57.502  32.161  1.00 30.33           O  
ATOM   1282  C3*   U A  59      71.352  58.357  34.366  1.00 30.46           C  
ATOM   1283  O3*   U A  59      70.546  58.855  35.497  1.00 31.17           O  
ATOM   1284  C2*   U A  59      71.629  56.872  34.435  1.00 34.51           C  
ATOM   1285  O2*   U A  59      70.529  56.200  35.034  1.00 33.82           O  
ATOM   1286  C1*   U A  59      71.808  56.483  32.934  1.00 34.05           C  
ATOM   1287  N1    U A  59      73.191  56.244  32.453  1.00 30.28           N  
ATOM   1288  C2    U A  59      73.828  55.055  32.883  1.00 33.87           C  
ATOM   1289  O2    U A  59      73.232  54.201  33.510  1.00 40.51           O  
ATOM   1290  N3    U A  59      75.160  54.920  32.534  1.00 31.37           N  
ATOM   1291  C4    U A  59      75.886  55.811  31.767  1.00 33.85           C  
ATOM   1292  O4    U A  59      76.986  55.468  31.336  1.00 34.46           O  
ATOM   1293  C5    U A  59      75.128  57.004  31.315  1.00 30.62           C  
ATOM   1294  C6    U A  59      73.842  57.148  31.666  1.00 24.94           C  
ATOM   1295  P     C A  60      71.267  59.539  36.803  1.00 33.56           P  
ATOM   1296  O1P   C A  60      70.297  59.837  37.820  1.00 32.77           O  
ATOM   1297  O2P   C A  60      72.180  60.638  36.294  1.00 34.80           O  
ATOM   1298  O5*   C A  60      72.351  58.459  37.303  1.00 29.47           O  
ATOM   1299  C5*   C A  60      71.983  57.223  37.962  1.00 31.76           C  
ATOM   1300  C4*   C A  60      73.252  56.436  38.285  1.00 32.79           C  
ATOM   1301  O4*   C A  60      74.024  56.360  37.069  1.00 34.84           O  
ATOM   1302  C3*   C A  60      74.192  57.080  39.322  1.00 35.07           C  
ATOM   1303  O3*   C A  60      74.832  56.043  40.028  1.00 38.94           O  
ATOM   1304  C2*   C A  60      75.220  57.789  38.455  1.00 34.97           C  
ATOM   1305  O2*   C A  60      76.489  58.039  39.019  1.00 36.28           O  
ATOM   1306  C1*   C A  60      75.336  56.786  37.306  1.00 32.64           C  
ATOM   1307  N1    C A  60      75.849  57.370  36.090  1.00 32.74           N  
ATOM   1308  C2    C A  60      76.958  56.775  35.475  1.00 33.30           C  
ATOM   1309  O2    C A  60      77.484  55.783  36.004  1.00 37.02           O  
ATOM   1310  N3    C A  60      77.445  57.308  34.368  1.00 27.75           N  
ATOM   1311  C4    C A  60      76.919  58.410  33.866  1.00 31.44           C  
ATOM   1312  N4    C A  60      77.491  58.943  32.854  1.00 31.16           N  
ATOM   1313  C5    C A  60      75.767  59.040  34.441  1.00 31.34           C  
ATOM   1314  C6    C A  60      75.274  58.489  35.549  1.00 34.59           C  
ATOM   1315  P     C A  61      74.243  55.584  41.447  1.00 42.61           P  
ATOM   1316  O1P   C A  61      75.151  54.490  41.830  1.00 44.00           O  
ATOM   1317  O2P   C A  61      72.781  55.398  41.440  1.00 37.71           O  
ATOM   1318  O5*   C A  61      74.477  56.853  42.390  1.00 39.76           O  
ATOM   1319  C5*   C A  61      75.717  57.059  43.037  1.00 43.07           C  
ATOM   1320  C4*   C A  61      75.643  58.316  43.840  1.00 39.59           C  
ATOM   1321  O4*   C A  61      75.807  59.456  42.959  1.00 37.14           O  
ATOM   1322  C3*   C A  61      74.297  58.531  44.524  1.00 38.59           C  
ATOM   1323  O3*   C A  61      74.209  57.799  45.784  1.00 39.28           O  
ATOM   1324  C2*   C A  61      74.292  60.040  44.702  1.00 37.97           C  
ATOM   1325  O2*   C A  61      75.174  60.397  45.747  1.00 38.59           O  
ATOM   1326  C1*   C A  61      74.974  60.515  43.421  1.00 36.12           C  
ATOM   1327  N1    C A  61      74.069  60.960  42.325  1.00 33.98           N  
ATOM   1328  C2    C A  61      73.361  62.180  42.471  1.00 33.04           C  
ATOM   1329  O2    C A  61      73.490  62.839  43.524  1.00 35.58           O  
ATOM   1330  N3    C A  61      72.551  62.600  41.455  1.00 35.45           N  
ATOM   1331  C4    C A  61      72.421  61.838  40.349  1.00 35.34           C  
ATOM   1332  N4    C A  61      71.550  62.237  39.373  1.00 32.83           N  
ATOM   1333  C5    C A  61      73.147  60.623  40.187  1.00 35.42           C  
ATOM   1334  C6    C A  61      73.938  60.221  41.188  1.00 31.03           C  
ATOM   1335  P     A A  62      72.777  57.319  46.334  1.00 42.33           P  
ATOM   1336  O1P   A A  62      73.023  56.390  47.455  1.00 41.97           O  
ATOM   1337  O2P   A A  62      71.960  56.820  45.161  1.00 42.66           O  
ATOM   1338  O5*   A A  62      72.104  58.670  46.874  1.00 43.04           O  
ATOM   1339  C5*   A A  62      72.719  59.460  47.936  1.00 39.06           C  
ATOM   1340  C4*   A A  62      71.882  60.684  48.244  1.00 39.69           C  
ATOM   1341  O4*   A A  62      71.968  61.660  47.165  1.00 39.39           O  
ATOM   1342  C3*   A A  62      70.391  60.442  48.405  1.00 40.80           C  
ATOM   1343  O3*   A A  62      70.103  59.945  49.700  1.00 43.31           O  
ATOM   1344  C2*   A A  62      69.808  61.827  48.120  1.00 41.10           C  
ATOM   1345  O2*   A A  62      70.016  62.793  49.130  1.00 43.43           O  
ATOM   1346  C1*   A A  62      70.711  62.312  46.995  1.00 38.85           C  
ATOM   1347  N9    A A  62      70.182  62.021  45.654  1.00 32.58           N  
ATOM   1348  C8    A A  62      70.560  61.004  44.827  1.00 33.66           C  
ATOM   1349  N7    A A  62      70.002  61.038  43.650  1.00 32.82           N  
ATOM   1350  C5    A A  62      69.158  62.140  43.707  1.00 31.78           C  
ATOM   1351  C6    A A  62      68.304  62.715  42.759  1.00 29.57           C  
ATOM   1352  N6    A A  62      68.170  62.255  41.478  1.00 29.13           N  
ATOM   1353  N1    A A  62      67.590  63.788  43.135  1.00 31.96           N  
ATOM   1354  C2    A A  62      67.754  64.268  44.373  1.00 33.29           C  
ATOM   1355  N3    A A  62      68.548  63.835  45.349  1.00 32.89           N  
ATOM   1356  C4    A A  62      69.239  62.747  44.945  1.00 35.51           C  
ATOM   1357  P     C A  63      68.913  58.881  49.918  1.00 43.35           P  
ATOM   1358  O1P   C A  63      69.009  58.532  51.351  1.00 49.74           O  
ATOM   1359  O2P   C A  63      68.970  57.790  48.905  1.00 46.77           O  
ATOM   1360  O5*   C A  63      67.599  59.762  49.769  1.00 45.50           O  
ATOM   1361  C5*   C A  63      67.371  60.836  50.703  1.00 44.28           C  
ATOM   1362  C4*   C A  63      66.263  61.750  50.244  1.00 44.77           C  
ATOM   1363  O4*   C A  63      66.669  62.512  49.079  1.00 41.48           O  
ATOM   1364  C3*   C A  63      64.958  61.089  49.850  1.00 46.16           C  
ATOM   1365  O3*   C A  63      64.160  60.771  51.005  1.00 51.01           O  
ATOM   1366  C2*   C A  63      64.325  62.154  48.972  1.00 43.55           C  
ATOM   1367  O2*   C A  63      63.764  63.209  49.763  1.00 48.54           O  
ATOM   1368  C1*   C A  63      65.541  62.713  48.233  1.00 39.95           C  
ATOM   1369  N1    C A  63      65.788  62.041  46.918  1.00 35.72           N  
ATOM   1370  C2    C A  63      65.155  62.553  45.793  1.00 32.68           C  
ATOM   1371  O2    C A  63      64.440  63.541  45.950  1.00 36.07           O  
ATOM   1372  N3    C A  63      65.328  61.946  44.571  1.00 34.10           N  
ATOM   1373  C4    C A  63      66.113  60.839  44.485  1.00 31.53           C  
ATOM   1374  N4    C A  63      66.253  60.229  43.282  1.00 34.07           N  
ATOM   1375  C5    C A  63      66.783  60.308  45.617  1.00 35.81           C  
ATOM   1376  C6    C A  63      66.603  60.941  46.809  1.00 36.37           C  
ATOM   1377  P     A A  64      63.129  59.542  50.943  1.00 57.75           P  
ATOM   1378  O1P   A A  64      62.528  59.402  52.300  1.00 58.22           O  
ATOM   1379  O2P   A A  64      63.807  58.359  50.338  1.00 54.03           O  
ATOM   1380  O5*   A A  64      62.008  60.100  49.971  1.00 54.48           O  
ATOM   1381  C5*   A A  64      61.194  61.198  50.398  1.00 56.10           C  
ATOM   1382  C4*   A A  64      60.204  61.569  49.339  1.00 54.58           C  
ATOM   1383  O4*   A A  64      60.917  62.255  48.276  1.00 53.92           O  
ATOM   1384  C3*   A A  64      59.459  60.433  48.632  1.00 56.48           C  
ATOM   1385  O3*   A A  64      58.374  59.780  49.384  1.00 56.71           O  
ATOM   1386  C2*   A A  64      59.036  61.152  47.349  1.00 53.00           C  
ATOM   1387  O2*   A A  64      57.994  62.076  47.565  1.00 56.42           O  
ATOM   1388  C1*   A A  64      60.283  61.983  47.026  1.00 51.87           C  
ATOM   1389  N9    A A  64      61.176  61.174  46.178  1.00 45.48           N  
ATOM   1390  C8    A A  64      62.220  60.374  46.567  1.00 42.54           C  
ATOM   1391  N7    A A  64      62.766  59.702  45.585  1.00 42.30           N  
ATOM   1392  C5    A A  64      62.051  60.122  44.466  1.00 40.22           C  
ATOM   1393  C6    A A  64      62.152  59.783  43.119  1.00 36.50           C  
ATOM   1394  N6    A A  64      63.082  58.963  42.657  1.00 34.41           N  
ATOM   1395  N1    A A  64      61.261  60.326  42.262  1.00 37.69           N  
ATOM   1396  C2    A A  64      60.346  61.178  42.745  1.00 36.62           C  
ATOM   1397  N3    A A  64      60.169  61.600  43.991  1.00 36.71           N  
ATOM   1398  C4    A A  64      61.070  61.018  44.815  1.00 39.74           C  
ATOM   1399  P     G A  65      58.068  58.186  49.152  1.00 58.45           P  
ATOM   1400  O1P   G A  65      57.104  57.682  50.155  1.00 62.97           O  
ATOM   1401  O2P   G A  65      59.328  57.423  48.941  1.00 61.74           O  
ATOM   1402  O5*   G A  65      57.302  58.169  47.766  1.00 57.06           O  
ATOM   1403  C5*   G A  65      56.297  59.121  47.518  1.00 49.45           C  
ATOM   1404  C4*   G A  65      55.988  59.170  46.060  1.00 45.91           C  
ATOM   1405  O4*   G A  65      57.113  59.757  45.358  1.00 43.81           O  
ATOM   1406  C3*   G A  65      55.711  57.868  45.316  1.00 45.00           C  
ATOM   1407  O3*   G A  65      54.366  57.370  45.509  1.00 47.48           O  
ATOM   1408  C2*   G A  65      55.920  58.322  43.886  1.00 41.79           C  
ATOM   1409  O2*   G A  65      54.776  59.084  43.491  1.00 39.93           O  
ATOM   1410  C1*   G A  65      57.142  59.262  44.032  1.00 42.23           C  
ATOM   1411  N9    G A  65      58.373  58.489  43.875  1.00 39.23           N  
ATOM   1412  C8    G A  65      59.245  58.086  44.861  1.00 37.65           C  
ATOM   1413  N7    G A  65      60.189  57.305  44.420  1.00 35.71           N  
ATOM   1414  C5    G A  65      59.942  57.213  43.050  1.00 38.74           C  
ATOM   1415  C6    G A  65      60.647  56.509  42.009  1.00 36.19           C  
ATOM   1416  O6    G A  65      61.690  55.776  42.101  1.00 38.10           O  
ATOM   1417  N1    G A  65      60.040  56.710  40.766  1.00 37.09           N  
ATOM   1418  C2    G A  65      58.920  57.482  40.550  1.00 37.20           C  
ATOM   1419  N2    G A  65      58.460  57.534  39.290  1.00 36.90           N  
ATOM   1420  N3    G A  65      58.288  58.152  41.498  1.00 33.44           N  
ATOM   1421  C4    G A  65      58.838  57.965  42.703  1.00 36.71           C  
ATOM   1422  P     A A  66      54.074  55.782  45.410  1.00 48.51           P  
ATOM   1423  O1P   A A  66      52.701  55.474  45.862  1.00 50.89           O  
ATOM   1424  O2P   A A  66      55.214  55.033  46.008  1.00 48.73           O  
ATOM   1425  O5*   A A  66      54.173  55.469  43.863  1.00 45.90           O  
ATOM   1426  C5*   A A  66      53.339  56.143  42.940  1.00 46.19           C  
ATOM   1427  C4*   A A  66      53.671  55.677  41.542  1.00 46.27           C  
ATOM   1428  O4*   A A  66      55.037  56.058  41.206  1.00 44.63           O  
ATOM   1429  C3*   A A  66      53.666  54.171  41.354  1.00 44.76           C  
ATOM   1430  O3*   A A  66      52.326  53.680  41.155  1.00 46.57           O  
ATOM   1431  C2*   A A  66      54.535  54.017  40.112  1.00 42.31           C  
ATOM   1432  O2*   A A  66      53.773  54.358  38.980  1.00 44.56           O  
ATOM   1433  C1*   A A  66      55.595  55.117  40.311  1.00 40.50           C  
ATOM   1434  N9    A A  66      56.838  54.618  40.902  1.00 38.52           N  
ATOM   1435  C8    A A  66      57.216  54.722  42.216  1.00 37.83           C  
ATOM   1436  N7    A A  66      58.365  54.141  42.488  1.00 38.64           N  
ATOM   1437  C5    A A  66      58.784  53.639  41.260  1.00 37.61           C  
ATOM   1438  C6    A A  66      59.943  52.898  40.877  1.00 36.94           C  
ATOM   1439  N6    A A  66      60.924  52.566  41.728  1.00 34.85           N  
ATOM   1440  N1    A A  66      60.042  52.511  39.596  1.00 36.33           N  
ATOM   1441  C2    A A  66      59.056  52.850  38.752  1.00 37.97           C  
ATOM   1442  N3    A A  66      57.922  53.552  38.986  1.00 38.42           N  
ATOM   1443  C4    A A  66      57.850  53.916  40.276  1.00 37.28           C  
ATOM   1444  P     A A  67      51.950  52.180  41.620  1.00 46.84           P  
ATOM   1445  O1P   A A  67      50.472  52.017  41.466  1.00 50.38           O  
ATOM   1446  O2P   A A  67      52.569  51.904  42.938  1.00 43.27           O  
ATOM   1447  O5*   A A  67      52.644  51.252  40.553  1.00 44.33           O  
ATOM   1448  C5*   A A  67      52.276  51.321  39.167  1.00 45.20           C  
ATOM   1449  C4*   A A  67      53.319  50.645  38.320  1.00 43.61           C  
ATOM   1450  O4*   A A  67      54.601  51.353  38.446  1.00 42.36           O  
ATOM   1451  C3*   A A  67      53.636  49.208  38.725  1.00 41.30           C  
ATOM   1452  O3*   A A  67      52.690  48.291  38.163  1.00 46.52           O  
ATOM   1453  C2*   A A  67      55.053  49.024  38.157  1.00 42.21           C  
ATOM   1454  O2*   A A  67      55.039  48.819  36.759  1.00 41.99           O  
ATOM   1455  C1*   A A  67      55.669  50.417  38.391  1.00 41.64           C  
ATOM   1456  N9    A A  67      56.441  50.432  39.643  1.00 39.97           N  
ATOM   1457  C8    A A  67      56.129  50.866  40.893  1.00 40.04           C  
ATOM   1458  N7    A A  67      57.062  50.618  41.788  1.00 37.46           N  
ATOM   1459  C5    A A  67      58.077  50.015  41.060  1.00 37.07           C  
ATOM   1460  C6    A A  67      59.348  49.481  41.435  1.00 37.10           C  
ATOM   1461  N6    A A  67      59.864  49.515  42.677  1.00 37.48           N  
ATOM   1462  N1    A A  67      60.076  48.895  40.483  1.00 37.12           N  
ATOM   1463  C2    A A  67      59.577  48.842  39.230  1.00 42.46           C  
ATOM   1464  N3    A A  67      58.415  49.313  38.758  1.00 39.78           N  
ATOM   1465  C4    A A  67      57.712  49.894  39.743  1.00 37.61           C  
ATOM   1466  P     U A  68      52.371  46.894  38.917  1.00 42.45           P  
ATOM   1467  O1P   U A  68      51.176  46.256  38.300  1.00 49.18           O  
ATOM   1468  O2P   U A  68      52.399  47.127  40.381  1.00 41.13           O  
ATOM   1469  O5*   U A  68      53.625  45.986  38.529  1.00 44.81           O  
ATOM   1470  C5*   U A  68      53.889  45.709  37.150  1.00 43.19           C  
ATOM   1471  C4*   U A  68      55.168  44.940  37.012  1.00 43.76           C  
ATOM   1472  O4*   U A  68      56.296  45.792  37.351  1.00 42.13           O  
ATOM   1473  C3*   U A  68      55.315  43.758  37.955  1.00 41.59           C  
ATOM   1474  O3*   U A  68      54.653  42.598  37.490  1.00 42.60           O  
ATOM   1475  C2*   U A  68      56.825  43.573  37.996  1.00 38.55           C  
ATOM   1476  O2*   U A  68      57.255  43.012  36.776  1.00 36.20           O  
ATOM   1477  C1*   U A  68      57.286  45.028  38.018  1.00 40.75           C  
ATOM   1478  N1    U A  68      57.451  45.543  39.385  1.00 40.94           N  
ATOM   1479  C2    U A  68      58.636  45.212  40.015  1.00 38.74           C  
ATOM   1480  O2    U A  68      59.470  44.485  39.490  1.00 38.70           O  
ATOM   1481  N3    U A  68      58.804  45.747  41.261  1.00 37.87           N  
ATOM   1482  C4    U A  68      57.924  46.564  41.939  1.00 41.13           C  
ATOM   1483  O4    U A  68      58.175  46.857  43.097  1.00 38.53           O  
ATOM   1484  C5    U A  68      56.693  46.846  41.241  1.00 42.38           C  
ATOM   1485  C6    U A  68      56.502  46.326  40.013  1.00 41.40           C  
ATOM   1486  P     U A  69      54.151  41.542  38.559  1.00 39.53           P  
ATOM   1487  O1P   U A  69      53.485  40.463  37.781  1.00 45.74           O  
ATOM   1488  O2P   U A  69      53.441  42.272  39.591  1.00 36.81           O  
ATOM   1489  O5*   U A  69      55.409  40.901  39.291  1.00 37.79           O  
ATOM   1490  C5*   U A  69      56.364  40.117  38.558  1.00 37.34           C  
ATOM   1491  C4*   U A  69      57.582  39.880  39.400  1.00 37.79           C  
ATOM   1492  O4*   U A  69      58.241  41.150  39.731  1.00 36.27           O  
ATOM   1493  C3*   U A  69      57.300  39.247  40.746  1.00 37.91           C  
ATOM   1494  O3*   U A  69      57.092  37.832  40.625  1.00 37.28           O  
ATOM   1495  C2*   U A  69      58.553  39.630  41.533  1.00 35.61           C  
ATOM   1496  O2*   U A  69      59.645  38.794  41.156  1.00 34.07           O  
ATOM   1497  C1*   U A  69      58.796  41.050  41.019  1.00 34.46           C  
ATOM   1498  N1    U A  69      58.186  42.082  41.876  1.00 33.94           N  
ATOM   1499  C2    U A  69      58.871  42.410  43.052  1.00 34.89           C  
ATOM   1500  O2    U A  69      59.947  41.899  43.373  1.00 36.74           O  
ATOM   1501  N3    U A  69      58.261  43.334  43.836  1.00 35.12           N  
ATOM   1502  C4    U A  69      57.071  43.985  43.603  1.00 37.89           C  
ATOM   1503  O4    U A  69      56.665  44.808  44.432  1.00 37.56           O  
ATOM   1504  C5    U A  69      56.419  43.627  42.355  1.00 38.56           C  
ATOM   1505  C6    U A  69      56.991  42.702  41.549  1.00 35.52           C  
ATOM   1506  P     C A  70      56.232  37.065  41.726  1.00 41.13           P  
ATOM   1507  O1P   C A  70      55.792  35.757  41.056  1.00 42.55           O  
ATOM   1508  O2P   C A  70      55.233  37.976  42.338  1.00 39.43           O  
ATOM   1509  O5*   C A  70      57.233  36.693  42.902  1.00 36.56           O  
ATOM   1510  C5*   C A  70      58.371  35.817  42.663  1.00 39.97           C  
ATOM   1511  C4*   C A  70      59.348  35.925  43.805  1.00 42.01           C  
ATOM   1512  O4*   C A  70      59.923  37.260  43.909  1.00 41.84           O  
ATOM   1513  C3*   C A  70      58.743  35.702  45.172  1.00 41.90           C  
ATOM   1514  O3*   C A  70      58.575  34.329  45.408  1.00 45.51           O  
ATOM   1515  C2*   C A  70      59.765  36.355  46.087  1.00 41.68           C  
ATOM   1516  O2*   C A  70      60.879  35.486  46.166  1.00 42.55           O  
ATOM   1517  C1*   C A  70      60.163  37.591  45.267  1.00 36.50           C  
ATOM   1518  N1    C A  70      59.360  38.780  45.597  1.00 39.49           N  
ATOM   1519  C2    C A  70      59.805  39.676  46.617  1.00 39.47           C  
ATOM   1520  O2    C A  70      60.879  39.451  47.194  1.00 41.25           O  
ATOM   1521  N3    C A  70      59.047  40.767  46.935  1.00 37.96           N  
ATOM   1522  C4    C A  70      57.901  40.982  46.285  1.00 40.08           C  
ATOM   1523  N4    C A  70      57.169  42.038  46.653  1.00 40.99           N  
ATOM   1524  C5    C A  70      57.445  40.109  45.229  1.00 39.06           C  
ATOM   1525  C6    C A  70      58.187  39.026  44.934  1.00 35.59           C  
ATOM   1526  P     G A  71      57.514  33.856  46.509  1.00 51.79           P  
ATOM   1527  O1P   G A  71      57.674  32.367  46.524  1.00 50.10           O  
ATOM   1528  O2P   G A  71      56.203  34.485  46.254  1.00 46.01           O  
ATOM   1529  O5*   G A  71      58.085  34.416  47.884  1.00 49.09           O  
ATOM   1530  C5*   G A  71      59.164  33.703  48.505  1.00 56.68           C  
ATOM   1531  C4*   G A  71      59.688  34.434  49.701  1.00 58.95           C  
ATOM   1532  O4*   G A  71      60.125  35.757  49.305  1.00 57.94           O  
ATOM   1533  C3*   G A  71      58.716  34.693  50.839  1.00 63.78           C  
ATOM   1534  O3*   G A  71      58.514  33.580  51.699  1.00 70.59           O  
ATOM   1535  C2*   G A  71      59.392  35.845  51.559  1.00 63.20           C  
ATOM   1536  O2*   G A  71      60.455  35.370  52.375  1.00 64.19           O  
ATOM   1537  C1*   G A  71      59.915  36.661  50.374  1.00 59.19           C  
ATOM   1538  N9    G A  71      58.897  37.620  49.946  1.00 58.54           N  
ATOM   1539  C8    G A  71      58.004  37.488  48.897  1.00 56.09           C  
ATOM   1540  N7    G A  71      57.207  38.519  48.769  1.00 54.09           N  
ATOM   1541  C5    G A  71      57.597  39.383  49.788  1.00 53.61           C  
ATOM   1542  C6    G A  71      57.114  40.666  50.148  1.00 52.69           C  
ATOM   1543  O6    G A  71      56.184  41.315  49.653  1.00 52.69           O  
ATOM   1544  N1    G A  71      57.827  41.196  51.223  1.00 52.34           N  
ATOM   1545  C2    G A  71      58.844  40.566  51.884  1.00 53.05           C  
ATOM   1546  N2    G A  71      59.399  41.251  52.885  1.00 55.23           N  
ATOM   1547  N3    G A  71      59.291  39.360  51.585  1.00 54.72           N  
ATOM   1548  C4    G A  71      58.636  38.835  50.526  1.00 55.38           C  
ATOM   1549  P     C A  72      57.094  33.402  52.426  1.00 73.67           P  
ATOM   1550  O1P   C A  72      57.177  32.132  53.187  1.00 76.79           O  
ATOM   1551  O2P   C A  72      56.011  33.588  51.427  1.00 75.72           O  
ATOM   1552  O5*   C A  72      57.017  34.599  53.469  1.00 74.00           O  
ATOM   1553  C5*   C A  72      57.817  34.563  54.660  1.00 75.03           C  
ATOM   1554  C4*   C A  72      57.653  35.834  55.450  1.00 74.38           C  
ATOM   1555  O4*   C A  72      57.964  36.953  54.582  1.00 74.64           O  
ATOM   1556  C3*   C A  72      56.261  36.155  55.979  1.00 76.09           C  
ATOM   1557  O3*   C A  72      55.932  35.502  57.209  1.00 78.85           O  
ATOM   1558  C2*   C A  72      56.302  37.671  56.139  1.00 74.55           C  
ATOM   1559  O2*   C A  72      56.872  38.146  57.338  1.00 75.06           O  
ATOM   1560  C1*   C A  72      57.184  38.082  54.963  1.00 72.55           C  
ATOM   1561  N1    C A  72      56.343  38.502  53.829  1.00 68.80           N  
ATOM   1562  C2    C A  72      55.796  39.801  53.838  1.00 66.61           C  
ATOM   1563  O2    C A  72      56.063  40.568  54.791  1.00 63.14           O  
ATOM   1564  N3    C A  72      54.993  40.180  52.815  1.00 64.67           N  
ATOM   1565  C4    C A  72      54.734  39.331  51.818  1.00 64.57           C  
ATOM   1566  N4    C A  72      53.935  39.749  50.832  1.00 63.54           N  
ATOM   1567  C5    C A  72      55.283  38.015  51.780  1.00 64.77           C  
ATOM   1568  C6    C A  72      56.077  37.647  52.792  1.00 67.37           C  
ATOM   1569  P     A A  73      54.383  35.178  57.550  1.00 79.81           P  
ATOM   1570  O1P   A A  73      54.333  34.483  58.872  1.00 81.24           O  
ATOM   1571  O2P   A A  73      53.785  34.517  56.355  1.00 79.38           O  
ATOM   1572  O5*   A A  73      53.739  36.617  57.766  1.00 77.14           O  
ATOM   1573  C5*   A A  73      54.204  37.443  58.849  1.00 75.83           C  
ATOM   1574  C4*   A A  73      53.508  38.772  58.837  1.00 74.32           C  
ATOM   1575  O4*   A A  73      53.912  39.514  57.654  1.00 72.81           O  
ATOM   1576  C3*   A A  73      51.987  38.727  58.750  1.00 73.77           C  
ATOM   1577  O3*   A A  73      51.332  38.492  59.999  1.00 76.93           O  
ATOM   1578  C2*   A A  73      51.667  40.101  58.178  1.00 71.67           C  
ATOM   1579  O2*   A A  73      51.704  41.127  59.141  1.00 70.42           O  
ATOM   1580  C1*   A A  73      52.823  40.304  57.196  1.00 69.17           C  
ATOM   1581  N9    A A  73      52.434  39.853  55.856  1.00 64.04           N  
ATOM   1582  C8    A A  73      52.792  38.695  55.206  1.00 60.25           C  
ATOM   1583  N7    A A  73      52.256  38.571  54.016  1.00 58.41           N  
ATOM   1584  C5    A A  73      51.493  39.729  53.869  1.00 57.85           C  
ATOM   1585  C6    A A  73      50.683  40.203  52.812  1.00 56.02           C  
ATOM   1586  N6    A A  73      50.516  39.541  51.654  1.00 55.49           N  
ATOM   1587  N1    A A  73      50.052  41.388  52.979  1.00 53.67           N  
ATOM   1588  C2    A A  73      50.233  42.052  54.136  1.00 55.48           C  
ATOM   1589  N3    A A  73      50.975  41.713  55.198  1.00 59.62           N  
ATOM   1590  C4    A A  73      51.589  40.525  54.993  1.00 59.80           C  
ATOM   1591  P     C A  74      49.934  37.680  60.030  1.00 78.48           P  
ATOM   1592  O1P   C A  74      49.485  37.580  61.446  1.00 80.11           O  
ATOM   1593  O2P   C A  74      50.090  36.441  59.222  1.00 77.06           O  
ATOM   1594  O5*   C A  74      48.914  38.636  59.275  1.00 77.54           O  
ATOM   1595  C5*   C A  74      48.529  39.898  59.843  1.00 79.25           C  
ATOM   1596  C4*   C A  74      47.463  40.541  58.985  1.00 80.49           C  
ATOM   1597  O4*   C A  74      48.030  40.879  57.690  1.00 79.75           O  
ATOM   1598  C3*   C A  74      46.291  39.627  58.660  1.00 81.62           C  
ATOM   1599  O3*   C A  74      45.292  39.640  59.669  1.00 82.15           O  
ATOM   1600  C2*   C A  74      45.778  40.187  57.340  1.00 80.42           C  
ATOM   1601  O2*   C A  74      44.932  41.302  57.547  1.00 82.05           O  
ATOM   1602  C1*   C A  74      47.080  40.621  56.662  1.00 79.08           C  
ATOM   1603  N1    C A  74      47.631  39.588  55.748  1.00 76.36           N  
ATOM   1604  C2    C A  74      47.204  39.557  54.400  1.00 74.77           C  
ATOM   1605  O2    C A  74      46.399  40.404  54.002  1.00 75.12           O  
ATOM   1606  N3    C A  74      47.694  38.606  53.571  1.00 73.87           N  
ATOM   1607  C4    C A  74      48.588  37.723  54.024  1.00 73.03           C  
ATOM   1608  N4    C A  74      49.055  36.813  53.167  1.00 72.80           N  
ATOM   1609  C5    C A  74      49.046  37.734  55.376  1.00 73.50           C  
ATOM   1610  C6    C A  74      48.543  38.671  56.196  1.00 74.97           C  
ATOM   1611  P     C A  75      44.492  38.295  59.994  1.00 83.00           P  
ATOM   1612  O1P   C A  75      43.435  38.637  60.974  1.00 84.71           O  
ATOM   1613  O2P   C A  75      45.478  37.229  60.303  1.00 81.99           O  
ATOM   1614  O5*   C A  75      43.802  37.890  58.618  1.00 84.55           O  
ATOM   1615  C5*   C A  75      42.993  38.825  57.874  1.00 85.33           C  
ATOM   1616  C4*   C A  75      42.564  38.194  56.568  1.00 86.51           C  
ATOM   1617  O4*   C A  75      43.667  38.149  55.630  1.00 84.68           O  
ATOM   1618  C3*   C A  75      42.136  36.747  56.733  1.00 88.27           C  
ATOM   1619  O3*   C A  75      40.778  36.695  57.104  1.00 95.04           O  
ATOM   1620  C2*   C A  75      42.402  36.136  55.366  1.00 85.58           C  
ATOM   1621  O2*   C A  75      41.337  36.317  54.459  1.00 84.94           O  
ATOM   1622  C1*   C A  75      43.631  36.926  54.915  1.00 82.57           C  
ATOM   1623  N1    C A  75      44.883  36.201  55.170  1.00 78.44           N  
ATOM   1624  C2    C A  75      45.365  35.327  54.193  1.00 75.98           C  
ATOM   1625  O2    C A  75      44.719  35.169  53.149  1.00 74.80           O  
ATOM   1626  N3    C A  75      46.512  34.673  54.406  1.00 74.76           N  
ATOM   1627  C4    C A  75      47.176  34.846  55.543  1.00 75.41           C  
ATOM   1628  N4    C A  75      48.311  34.161  55.705  1.00 75.90           N  
ATOM   1629  C5    C A  75      46.709  35.719  56.563  1.00 75.64           C  
ATOM   1630  C6    C A  75      45.570  36.373  56.337  1.00 77.27           C  
ATOM   1631  P     A A  76      40.334  35.770  58.330  1.00 99.54           P  
ATOM   1632  O1P   A A  76      41.267  36.020  59.481  1.00 99.74           O  
ATOM   1633  O2P   A A  76      40.226  34.405  57.758  1.00 99.88           O  
ATOM   1634  O5*   A A  76      38.872  36.304  58.697  1.00100.19           O  
ATOM   1635  C5*   A A  76      38.666  37.596  59.323  1.00100.18           C  
ATOM   1636  C4*   A A  76      37.607  38.379  58.569  1.00100.19           C  
ATOM   1637  O4*   A A  76      36.479  37.500  58.278  1.00100.19           O  
ATOM   1638  C3*   A A  76      37.025  39.585  59.305  1.00100.19           C  
ATOM   1639  O3*   A A  76      36.428  40.400  58.274  1.00100.19           O  
ATOM   1640  C2*   A A  76      35.785  38.991  59.981  1.00100.19           C  
ATOM   1641  O2*   A A  76      34.780  39.945  60.283  1.00100.19           O  
ATOM   1642  C1*   A A  76      35.314  37.959  58.950  1.00100.19           C  
ATOM   1643  N9    A A  76      34.598  36.785  59.488  1.00100.19           N  
ATOM   1644  C8    A A  76      34.399  35.586  58.829  1.00100.19           C  
ATOM   1645  N7    A A  76      33.715  34.701  59.522  1.00100.19           N  
ATOM   1646  C5    A A  76      33.440  35.352  60.719  1.00100.19           C  
ATOM   1647  C6    A A  76      32.739  34.948  61.881  1.00100.19           C  
ATOM   1648  N6    A A  76      32.161  33.747  62.024  1.00100.19           N  
ATOM   1649  N1    A A  76      32.652  35.835  62.902  1.00100.19           N  
ATOM   1650  C2    A A  76      33.230  37.042  62.758  1.00100.19           C  
ATOM   1651  N3    A A  76      33.911  37.538  61.720  1.00100.19           N  
ATOM   1652  C4    A A  76      33.982  36.637  60.719  1.00100.19           C  
TER    1653        A A  76                                                      
HETATM 1654 MG    MG   590      70.566  35.530   1.665  1.00 65.44          MG  
HETATM 1655 MN   O4M   530      80.714  57.068  31.271  1.00 34.23          MN  
HETATM 1656  O4  O4M   530      80.684  55.382  30.181  1.00 37.94           O  
HETATM 1657  O3  O4M   530      81.503  58.060  29.728  1.00 34.54           O  
HETATM 1658  O2  O4M   530      82.554  56.601  31.916  1.00 40.13           O  
HETATM 1659  O1  O4M   530      78.882  57.524  30.607  1.00 37.65           O  
HETATM 1660 MG   MO5   510      57.346  47.575  47.279  1.00 62.77          MG  
HETATM 1661  O1  MO5   510      59.121  46.716  46.939  1.00 61.12           O  
HETATM 1662  O2  MO5   510      55.569  48.426  47.617  1.00 62.44           O  
HETATM 1663  O3  MO5   510      57.139  47.928  45.318  1.00 62.18           O  
HETATM 1664  O4  MO5   510      57.549  47.218  49.240  1.00 63.71           O  
HETATM 1665  O5  MO5   510      56.451  45.796  47.050  1.00 64.69           O  
HETATM 1666 MN   MN5   520      49.923  44.427  50.131  1.00 89.20          MN  
HETATM 1667  O1  MN5   520      50.566  42.548  49.894  1.00 90.06           O  
HETATM 1668  O2  MN5   520      49.279  46.306  50.367  1.00 88.19           O  
HETATM 1669  O3  MN5   520      48.828  43.854  51.702  1.00 89.34           O  
HETATM 1670  O4  MN5   520      51.020  45.000  48.559  1.00 89.35           O  
HETATM 1671  O5  MN5   520      48.379  44.049  48.917  1.00 88.93           O  
HETATM 1672 MG   MO3   540      77.110  64.307  25.357  1.00 40.08          MG  
HETATM 1673  O1  MO3   540      75.884  64.949  23.911  1.00 42.31           O  
HETATM 1674  O2  MO3   540      75.591  63.324  26.209  1.00 42.77           O  
HETATM 1675  O3  MO3   540      78.628  65.295  24.504  1.00 43.67           O  
HETATM 1676 MG   MO6   560      62.649  46.629  27.595  1.00 47.02          MG  
HETATM 1677  O1  MO6   560      63.354  44.755  27.511  1.00 43.55           O  
HETATM 1678  O2  MO6   560      61.942  48.502  27.674  1.00 41.87           O  
HETATM 1679  O3  MO6   560      62.352  46.432  29.566  1.00 44.71           O  
HETATM 1680  O4  MO6   560      62.949  46.821  25.628  1.00 43.68           O  
HETATM 1681  O5  MO6   560      64.498  47.320  27.946  1.00 41.71           O  
HETATM 1682  O6  MO6   560      60.803  45.945  27.247  1.00 46.26           O  
HETATM 1683 MG   MO6   570      73.331  43.321  11.207  1.00 50.39          MG  
HETATM 1684  O1  MO6   570      72.795  42.398   9.514  1.00 50.20           O  
HETATM 1685  O2  MO6   570      73.865  44.246  12.908  1.00 49.37           O  
HETATM 1686  O3  MO6   570      74.746  41.940  11.519  1.00 51.45           O  
HETATM 1687  O4  MO6   570      71.918  44.704  10.896  1.00 48.14           O  
HETATM 1688  O5  MO6   570      72.020  42.211  12.224  1.00 47.02           O  
HETATM 1689  O6  MO6   570      74.644  44.433  10.185  1.00 47.69           O  
HETATM 1690 MG   MO1   580      69.222  44.815  33.339  1.00 61.74          MG  
HETATM 1691  O1  MO1   580      68.372  45.072  31.544  1.00 46.52           O  
HETATM 1692 MN   MN5   550      72.301  48.513  33.894  1.00 56.51          MN  
HETATM 1693  O1  MN5   550      70.576  49.247  33.191  1.00 58.54           O  
HETATM 1694  O2  MN5   550      74.024  47.784  34.605  1.00 58.71           O  
HETATM 1695  O3  MN5   550      71.834  46.713  33.168  1.00 61.51           O  
HETATM 1696  O4  MN5   550      72.774  50.321  34.622  1.00 60.62           O  
HETATM 1697  O5  MN5   550      71.401  48.057  35.619  1.00 60.95           O  
HETATM 1698  O   HOH   101      65.235  47.736  24.306  1.00 29.97           O  
HETATM 1699  O   HOH   102      74.678  53.324  26.387  1.00 29.43           O  
HETATM 1700  O   HOH   103      79.647  66.543  30.502  1.00 36.21           O  
HETATM 1701  O   HOH   104      69.474  53.115  32.762  1.00 35.08           O  
HETATM 1702  O   HOH   105      77.803  59.348  29.070  1.00 34.91           O  
HETATM 1703  O   HOH   106      86.312  62.508  34.397  1.00 36.87           O  
HETATM 1704  O   HOH   107      69.798  47.380  19.420  1.00 35.43           O  
HETATM 1705  O   HOH   108      77.715  51.787  26.254  1.00 25.88           O  
HETATM 1706  O   HOH   109      66.697  54.043  19.991  1.00 38.52           O  
HETATM 1707  O   HOH   110      73.012  72.799  41.306  1.00 35.46           O  
HETATM 1708  O   HOH   111      84.966  51.999  36.727  1.00 48.17           O  
HETATM 1709  O   HOH   112      75.699  47.326  15.656  1.00 43.55           O  
HETATM 1710  O   HOH   113      61.911  39.313  42.834  1.00 38.72           O  
HETATM 1711  O   HOH   114      72.538  65.567  46.811  1.00 39.71           O  
HETATM 1712  O   HOH   115      64.957  57.362  35.588  1.00 36.28           O  
HETATM 1713  O   HOH   116      88.913  61.884  36.666  1.00 33.38           O  
HETATM 1714  O   HOH   117      77.430  47.049  24.346  1.00 38.25           O  
HETATM 1715  O   HOH   118      85.080  55.471  30.753  1.00 42.06           O  
HETATM 1716  O   HOH   119      73.126  42.774  26.088  1.00 39.89           O  
HETATM 1717  O   HOH   120      79.541  54.130  20.639  1.00 41.26           O  
HETATM 1718  O   HOH   121      75.971  56.472  28.176  1.00 40.47           O  
HETATM 1719  O   HOH   123      78.750  74.856  43.549  1.00 41.10           O  
HETATM 1720  O   HOH   124      59.778  66.139  42.613  1.00 45.92           O  
HETATM 1721  O   HOH   125      72.198  49.756  14.319  1.00 52.62           O  
HETATM 1722  O   HOH   126      68.821  65.008  47.959  1.00 45.12           O  
HETATM 1723  O   HOH   127      67.849  57.944  43.472  1.00 49.89           O  
HETATM 1724  O   HOH   128      67.631  51.090  33.685  1.00 39.79           O  
HETATM 1725  O   HOH   129      72.804  72.294  34.649  1.00 52.66           O  
HETATM 1726  O   HOH   130      71.861  49.100  17.322  1.00 38.23           O  
HETATM 1727  O   HOH   131      65.908  53.779  40.692  1.00 47.68           O  
HETATM 1728  O   HOH   132      54.150  42.092  45.268  1.00 43.73           O  
HETATM 1729  O   HOH   133      89.825  58.582  38.826  1.00 50.13           O  
HETATM 1730  O   HOH   134      84.722  60.865  36.177  1.00 45.21           O  
HETATM 1731  O   HOH   135      76.181  73.948  39.594  1.00 36.60           O  
HETATM 1732  O   HOH   136      77.944  64.645  31.617  1.00 50.38           O  
HETATM 1733  O   HOH   137      63.795  30.172  -8.994  1.00 51.48           O  
HETATM 1734  O   HOH   138      79.187  54.860  37.993  1.00 53.51           O  
HETATM 1735  O   HOH   139      65.438  43.393  28.966  1.00 49.20           O  
HETATM 1736  O   HOH   140      76.458  61.362  31.941  1.00 41.25           O  
HETATM 1737  O   HOH   141      65.955  45.704  31.155  1.00 37.31           O  
HETATM 1738  O   HOH   142      76.497  48.574  12.986  1.00 52.62           O  
HETATM 1739  O   HOH   143      77.696  63.153  34.024  1.00 41.50           O  
HETATM 1740  O   HOH   144      83.868  72.752  40.256  1.00 49.29           O  
HETATM 1741  O   HOH   145      83.766  59.152  28.946  1.00 39.50           O  
HETATM 1742  O   HOH   146      80.216  39.612  12.437  1.00 58.96           O  
HETATM 1743  O   HOH   147      74.386  47.570  12.387  1.00 53.73           O  
HETATM 1744  O   HOH   148      76.514  45.163  13.466  1.00 52.76           O  
HETATM 1745  O   HOH   149      63.032  41.158  41.088  1.00 48.34           O  
HETATM 1746  O   HOH   150      71.118  68.217  36.794  1.00 51.89           O  
HETATM 1747  O   HOH   152      81.091  65.399  26.812  1.00 44.48           O  
HETATM 1748  O   HOH   153      77.519  52.631  29.605  1.00 48.97           O  
HETATM 1749  O   HOH   155      59.968  42.360  37.108  1.00 49.01           O  
HETATM 1750  O   HOH   156      65.613  68.667  37.929  1.00 45.40           O  
HETATM 1751  O   HOH   157      56.060  63.049  38.179  1.00 49.08           O  
HETATM 1752  O   HOH   158      65.177  33.333   7.613  1.00 55.31           O  
HETATM 1753  O   HOH   159      75.172  67.087  29.927  1.00 54.05           O  
HETATM 1754  O   HOH   160      84.839  64.460  36.543  1.00 47.37           O  
HETATM 1755  O   HOH   161      84.226  60.935  26.655  1.00 51.90           O  
HETATM 1756  O   HOH   162      65.528  66.128  32.712  1.00 58.08           O  
HETATM 1757  O   HOH   163      54.801  39.353  47.303  1.00 50.47           O  
HETATM 1758  O   HOH   164      82.124  69.432  40.965  1.00 44.14           O  
HETATM 1759  O   HOH   166      76.820  63.025  29.658  1.00 45.12           O  
HETATM 1760  O   HOH   167      66.938  41.242  10.226  1.00 51.21           O  
HETATM 1761  O   HOH   168      85.023  60.379  30.965  1.00 38.68           O  
HETATM 1762  O   HOH   170      81.568  67.018  41.474  1.00 54.90           O  
HETATM 1763  O   HOH   171      83.672  64.839  39.206  1.00 48.85           O  
HETATM 1764  O   HOH   172      69.358  59.813  40.447  1.00 56.32           O  
HETATM 1765  O   HOH   173      65.402  39.517  46.250  1.00 42.43           O  
HETATM 1766  O   HOH   174      78.025  56.354  40.203  1.00 44.75           O  
HETATM 1767  O   HOH   175      72.640  54.075  19.103  1.00 56.38           O  
HETATM 1768  O   HOH   178      62.561  49.577   9.243  1.00 59.11           O  
HETATM 1769  O   HOH   181      67.851  57.255  35.234  1.00 39.17           O  
HETATM 1770  O   HOH   183      65.609  57.383  39.390  1.00 56.16           O  
HETATM 1771  O   HOH   184      77.652  63.370  43.463  1.00 46.31           O  
HETATM 1772  O   HOH   185      56.156  59.761  40.881  1.00 54.01           O  
HETATM 1773  O   HOH   186      68.030  57.904  37.829  1.00 54.73           O  
HETATM 1774  O   HOH   189      64.948  50.726  18.484  1.00 45.40           O  
HETATM 1775  O   HOH   191      58.581  70.881  37.367  1.00 43.35           O  
HETATM 1776  O   HOH   195      69.329  73.815  40.546  1.00 45.81           O  
HETATM 1777  O   HOH   196      71.092  44.843  13.706  1.00 46.18           O  
HETATM 1778  O   HOH   197      63.773  67.421  46.112  1.00 52.13           O  
HETATM 1779  O   HOH   200      79.526  44.399   2.237  1.00 56.28           O  
HETATM 1780  O   HOH   204      61.159  66.141  44.876  1.00 49.57           O  
HETATM 1781  O   HOH   205      55.921  58.490  38.100  1.00 52.97           O  
HETATM 1782  O   HOH   206      61.370  44.287  30.748  1.00 53.06           O  
HETATM 1783  O   HOH   208      72.463  66.452  30.831  1.00 58.61           O  
HETATM 1784  O   HOH   210      60.953  51.071  33.259  1.00 44.43           O  
HETATM 1785  O   HOH   214      55.561  30.912  50.683  1.00 57.65           O  
HETATM 1786  O   HOH   219      72.422  43.667  15.579  1.00 52.11           O  
HETATM 1787  O   HOH   222      65.477  55.377  26.488  1.00 40.16           O  
HETATM 1788  O   HOH   223      62.090  56.194  45.841  1.00 48.80           O  
HETATM 1789  O   HOH   226      60.948  49.649  25.176  1.00 34.03           O  
HETATM 1790  O   HOH   228      69.381  66.644  35.818  1.00 46.90           O  
HETATM 1791  O   HOH   230      66.314  66.119  35.359  1.00 45.34           O  
HETATM 1792  O   HOH   231      69.248  64.196  37.425  1.00 39.82           O  
HETATM 1793  O   HOH   233      67.490  67.214  38.107  1.00 42.96           O  
HETATM 1794  O   HOH   591      70.315  63.256  32.948  1.00 48.78           O  
HETATM 1795  O   HOH   592      70.278  67.253  46.788  1.00 49.58           O  
HETATM 1796  O   HOH   593      59.221  51.476  25.337  1.00 53.13           O  
HETATM 1797  O   HOH   596      36.759  35.239  56.874  1.00 54.15           O  
HETATM 1798  O   HOH   598      72.226  76.169  43.785  1.00 56.77           O  
HETATM 1799  O   HOH   602      79.271  66.367  28.116  1.00 57.39           O  
HETATM 1800  O   HOH   603      68.045  43.077   4.648  1.00 58.94           O  
HETATM 1801  O   HOH   608      52.188  37.918  39.395  1.00 56.79           O  
HETATM 1802  O   HOH   610      53.895  62.874  28.494  1.00 55.33           O  
HETATM 1803  O   HOH   611      70.166  44.500  35.929  1.00 59.11           O  
HETATM 1804  O   HOH   612      65.815  56.418  41.681  1.00 63.02           O  
HETATM 1805  O   HOH   613      84.445  61.099  24.272  1.00 53.13           O  
HETATM 1806  O   HOH   616      62.869  55.003  26.332  1.00 48.50           O  
HETATM 1807  O   HOH   618      50.840  53.124  44.965  1.00 59.63           O  
HETATM 1808  O   HOH   626      59.036  51.116  45.301  1.00 57.70           O  
HETATM 1809  O   HOH   627      51.263  39.619  37.323  1.00 55.93           O  
HETATM 1810  O   HOH   633      59.499  57.746  20.958  1.00 59.69           O  
HETATM 1811  O   HOH   635      77.328  58.675  26.658  1.00 51.93           O  
HETATM 1812  O   HOH   644      72.884  27.349  -8.052  1.00 57.57           O  
HETATM 1813  O   HOH   648      63.777  65.450  36.414  1.00 42.27           O  
HETATM 1814  O   HOH   657      72.947  30.361 -12.545  1.00 50.16           O  
HETATM 1815  O   HOH   662      57.486  68.555  38.998  1.00 41.15           O  
HETATM 1816  O   HOH   670      72.917  73.923  36.977  1.00 53.90           O  
HETATM 1817  O   HOH   671      82.577  50.441  36.557  1.00 57.25           O  
HETATM 1818  O   HOH   675      68.361  55.414  38.613  1.00 58.70           O  
HETATM 1819  O   HOH   688      64.284  30.562  -4.413  1.00 56.07           O  
HETATM 1820  O   HOH   690      69.590  42.478  12.375  1.00 58.33           O  
HETATM 1821  O   HOH   693      83.889  61.714  39.991  1.00 53.01           O  
HETATM 1822  O   HOH   707      61.422  49.192  45.932  1.00 56.92           O  
CONECT  181  195                                                                
CONECT  195  181  196  197  198                                                 
CONECT  196  195                                                                
CONECT  197  195                                                                
CONECT  198  195  199                                                           
CONECT  199  198  200                                                           
CONECT  200  199  201  202                                                      
CONECT  201  200  206                                                           
CONECT  202  200  203  204                                                      
CONECT  203  202  219                                                           
CONECT  204  202  205  206                                                      
CONECT  205  204                                                                
CONECT  206  201  204  207                                                      
CONECT  207  206  208  218                                                      
CONECT  208  207  209                                                           
CONECT  209  208  210                                                           
CONECT  210  209  211  218                                                      
CONECT  211  210  212  213                                                      
CONECT  212  211                                                                
CONECT  213  211  214                                                           
CONECT  214  213  215  217                                                      
CONECT  215  214  216                                                           
CONECT  216  215                                                                
CONECT  217  214  218                                                           
CONECT  218  207  210  217                                                      
CONECT  219  203                                                                
CONECT  309  324                                                                
CONECT  324  309  325  326  327                                                 
CONECT  325  324                                                                
CONECT  326  324                                                                
CONECT  327  324  328                                                           
CONECT  328  327  329                                                           
CONECT  329  328  330  331                                                      
CONECT  330  329  333                                                           
CONECT  331  329  332  334                                                      
CONECT  332  331  344                                                           
CONECT  333  330  334  336                                                      
CONECT  334  331  333  335                                                      
CONECT  335  334                                                                
CONECT  336  333  337  343                                                      
CONECT  337  336  338  339                                                      
CONECT  338  337                                                                
CONECT  339  337  340                                                           
CONECT  340  339  341  342                                                      
CONECT  341  340                                                                
CONECT  342  340  343                                                           
CONECT  343  336  342                                                           
CONECT  344  332  345  346  347                                                 
CONECT  345  344                                                                
CONECT  346  344                                                                
CONECT  347  344  348                                                           
CONECT  348  347  349                                                           
CONECT  349  348  350  351                                                      
CONECT  350  349  353                                                           
CONECT  351  349  352  354                                                      
CONECT  352  351  364                                                           
CONECT  353  350  354  356                                                      
CONECT  354  351  353  355                                                      
CONECT  355  354                                                                
CONECT  356  353  357  363                                                      
CONECT  357  356  358  359                                                      
CONECT  358  357                                                                
CONECT  359  357  360                                                           
CONECT  360  359  361  362                                                      
CONECT  361  360                                                                
CONECT  362  360  363                                                           
CONECT  363  356  362                                                           
CONECT  364  352                                                                
CONECT  531  543                                                                
CONECT  543  531  544  545  546                                                 
CONECT  544  543                                                                
CONECT  545  543                                                                
CONECT  546  543  547                                                           
CONECT  547  546  548                                                           
CONECT  548  547  549  550                                                      
CONECT  549  548  554                                                           
CONECT  550  548  551  552                                                      
CONECT  551  550  568                                                           
CONECT  552  550  553  554                                                      
CONECT  553  552                                                                
CONECT  554  549  552  555                                                      
CONECT  555  554  556  565                                                      
CONECT  556  555  557                                                           
CONECT  557  556  558                                                           
CONECT  558  557  559  565                                                      
CONECT  559  558  560  561                                                      
CONECT  560  559                                                                
CONECT  561  559  562                                                           
CONECT  562  561  563  564                                                      
CONECT  563  562  566  567                                                      
CONECT  564  562  565                                                           
CONECT  565  555  558  564                                                      
CONECT  566  563                                                                
CONECT  567  563                                                                
CONECT  568  551                                                                
CONECT  661  693                                                                
CONECT  675  676  680  683                                                      
CONECT  676  675  677  681                                                      
CONECT  677  676  678                                                           
CONECT  678  677  679  682                                                      
CONECT  679  678  680                                                           
CONECT  680  675  679                                                           
CONECT  681  676                                                                
CONECT  682  678                                                                
CONECT  683  675  684  689                                                      
CONECT  684  683  685  687                                                      
CONECT  685  684  686                                                           
CONECT  686  685                                                                
CONECT  687  684  688  690                                                      
CONECT  688  687  689  691                                                      
CONECT  689  683  688                                                           
CONECT  690  687  696                                                           
CONECT  691  688  692                                                           
CONECT  692  691  693                                                           
CONECT  693  661  692  694  695                                                 
CONECT  694  693                                                                
CONECT  695  693                                                                
CONECT  696  690                                                                
CONECT  704  716                                                                
CONECT  716  704  717  718  719                                                 
CONECT  717  716                                                                
CONECT  718  716                                                                
CONECT  719  716  720                                                           
CONECT  720  719  721                                                           
CONECT  721  720  722  723                                                      
CONECT  722  721  728                                                           
CONECT  723  721  724  725                                                      
CONECT  724  723  740                                                           
CONECT  725  723  726  728                                                      
CONECT  726  725  727                                                           
CONECT  727  726                                                                
CONECT  728  722  725  729                                                      
CONECT  729  728  730  739                                                      
CONECT  730  729  731                                                           
CONECT  731  730  732                                                           
CONECT  732  731  733  739                                                      
CONECT  733  732  734  735                                                      
CONECT  734  733                                                                
CONECT  735  733  736                                                           
CONECT  736  735  737  738                                                      
CONECT  737  736                                                                
CONECT  738  736  739                                                           
CONECT  739  729  732  738                                                      
CONECT  740  724                                                                
CONECT  770  820                                                                
CONECT  784  786  791  798                                                      
CONECT  785  786  797                                                           
CONECT  786  784  785  787                                                      
CONECT  787  786  788  789                                                      
CONECT  788  787                                                                
CONECT  789  787  790  795                                                      
CONECT  790  789  791  793                                                      
CONECT  791  784  790  792                                                      
CONECT  792  791                                                                
CONECT  793  790  794                                                           
CONECT  794  793  795                                                           
CONECT  795  789  794  811                                                      
CONECT  796  797                                                                
CONECT  797  785  796  798                                                      
CONECT  798  784  797  799                                                      
CONECT  799  798  800                                                           
CONECT  800  799  801                                                           
CONECT  801  800  802  806                                                      
CONECT  802  801  803  804                                                      
CONECT  803  802                                                                
CONECT  804  802  805                                                           
CONECT  805  804                                                                
CONECT  806  801  807                                                           
CONECT  807  806  808  809                                                      
CONECT  808  807                                                                
CONECT  809  807  810                                                           
CONECT  810  809                                                                
CONECT  811  795  812  817                                                      
CONECT  812  811  813  814                                                      
CONECT  813  812                                                                
CONECT  814  812  815  816                                                      
CONECT  815  814  823                                                           
CONECT  816  814  817  818                                                      
""" 

fake_nmr_file="""HEADER    RIBONUCLEIC ACID                        04-AUG-98   17RA              
TITLE     BRANCHPOINT HELIX FROM YEAST AND BINDING SITE FOR PHAGE               
TITLE    2 GA/MS2 COAT PROTEINS, NMR, 12 STRUCTURES                             
COMPND    MOL_ID: 1;                                                            
COMPND   2 MOLECULE: RNA;                                                       
COMPND   3 CHAIN: NULL;                                                         
COMPND   4 FRAGMENT: RBS AND START SITE FOR PHAGE GA REPLICASE GENE;            
COMPND   5 ENGINEERED: YES;                                                     
COMPND   6 MUTATION: A5U, A6U                                                   
SOURCE    MOL_ID: 1;                                                            
SOURCE   2 SYNTHETIC: YES;                                                      
SOURCE   3 ORGANISM_SCIENTIFIC: SACCHAROMYCES CEREVISIAE;                       
SOURCE   4 ORGANISM_COMMON: BAKER'S YEAST;                                      
SOURCE   5 CELLULAR_LOCATION: NUCLEUS;                                          
SOURCE   6 GENE: REPLICASE;                                                     
SOURCE   7 OTHER_DETAILS: IN VITRO SYNTHESIS FROM DNA TEMPLATE USING            
SOURCE   8 T7 RNA POLYMERASE. HAIRPIN CORRESPONDS TO NT -16 - +5 OF             
SOURCE   9 PHAGE GA REPLICASE AND THE YEAST PRE-MRNA BRANCHPOINT HELIX          
KEYWDS    BRANCHPOINT HELIX, PHAGE MS2, BULGE, BASE TRIPLE,                     
KEYWDS   2 RIBONUCLEIC ACID                                                     
EXPDTA    NMR, 12 STRUCTURES                                                    
AUTHOR    E.P.NIKONOWICZ,J.S.SMITH                                              
REVDAT   1   20-APR-99 17RA    0                                                
JRNL        AUTH   J.S.SMITH,E.P.NIKONOWICZ                                     
JRNL        TITL   NMR STRUCTURE AND DYNAMICS OF AN RNA MOTIF COMMON            
JRNL        TITL 2 TO THE SPLICEOSOME BRANCH-POINT HELIX AND THE                
JRNL        TITL 3 RNA-BINDING SITE FOR PHAGE GA COAT PROTEIN                   
JRNL        REF    BIOCHEMISTRY                  V.  37 13486 1998              
JRNL        REFN   ASTM BICHAW  US ISSN 0006-2960                 0033          
REMARK   1                                                                      
REMARK   2                                                                      
REMARK   2 RESOLUTION. NOT APPLICABLE.                                          
REMARK   3                                                                      
REMARK   3 REFINEMENT.                                                          
REMARK   3   PROGRAM     : X-PLOR 3.1                                           
REMARK   3   AUTHORS     : BRUNGER                                              
REMARK   3                                                                      
REMARK   3  OTHER REFINEMENT REMARKS: REFINEMENT DETAILS CAN BE FOUND           
REMARK   3  IN THE JRNL CITATION ABOVE.                                         
REMARK   4                                                                      
REMARK   4 17RA COMPLIES WITH FORMAT V. 2.3, 09-JULY-1998                       
REMARK   6                                                                      
REMARK   6 CALCULATIONS PERFORMED WITHOUT PHOSPHATE ON 5' TERMINAL              
REMARK   6 NUCLEOTIDE.                                                          
REMARK 210                                                                      
REMARK 210 EXPERIMENTAL DETAILS                                                 
REMARK 210  EXPERIMENT TYPE                : NMR                                
REMARK 210  TEMPERATURE           (KELVIN) : 298                                
REMARK 210  PH                             : 6.8                                
REMARK 210  IONIC STRENGTH                 : 20 MM                              
REMARK 210  PRESSURE                       : 1 ATM                              
REMARK 210  SAMPLE CONTENTS                : 10% H2O/90% D2O, 100% D2O          
REMARK 210                                                                      
REMARK 210  NMR EXPERIMENTS CONDUCTED      : 2D/3D NOESY, COSY, HETCOR          
REMARK 210  SPECTROMETER FIELD STRENGTH    : 500 MHZ                            
REMARK 210  SPECTROMETER MODEL             : AMX500                             
REMARK 210  SPECTROMETER MANUFACTURER      : BRUKER                             
REMARK 210                                                                      
REMARK 210  STRUCTURE DETERMINATION.                                            
REMARK 210   SOFTWARE USED                 : FELIX 950                          
REMARK 210   METHOD USED                   : SIMULATED ANNEALING                
REMARK 210                                                                      
REMARK 210 CONFORMERS, NUMBER CALCULATED   : 75                                 
REMARK 210 CONFORMERS, NUMBER SUBMITTED    : 12                                 
REMARK 210 CONFORMERS, SELECTION CRITERIA  : LEAST RESTRAINT VIOLATION          
REMARK 210                                                                      
REMARK 210 BEST REPRESENTATIVE CONFORMER IN THIS ENSEMBLE : NULL                
REMARK 210                                                                      
REMARK 210 REMARK: FURTHER DETAILS OF DATA ACQUISITION AND STRUCTURE            
REMARK 210 CALCULATION CAN BE FOUND IN THE METHODS SECTION OF THE               
REMARK 210 MANUSCRIPT LISTED ABOVE. MODELS 1-6 ARE MINIMIZED                    
REMARK 210 STRUCTURES CALCULATED WITHOUT BASE PAIR CONSTRAINTS FOR              
REMARK 210 A6, A7, AND U16. MODELS 7-12 ARE MINIMIZED STRUCTURES                
REMARK 210 CALCULATED USING A6-U16 BASE PAIR CONSTRAINTS AND AN A7              
REMARK 210 N1H-U16 O2 HYDROGEN BOND.                                            
REMARK 215                                                                      
REMARK 215 NMR STUDY                                                            
REMARK 215 THE COORDINATES IN THIS ENTRY WERE GENERATED FROM SOLUTION           
REMARK 215 NMR DATA.  PROTEIN DATA BANK CONVENTIONS REQUIRE THAT                
REMARK 215 CRYST1 AND SCALE RECORDS BE INCLUDED, BUT THE VALUES ON              
REMARK 215 THESE RECORDS ARE MEANINGLESS.                                       
REMARK 800                                                                      
REMARK 800 SITE                                                                 
REMARK 800 SITE_IDENTIFIER: AAU                                                 
REMARK 800 SITE_DESCRIPTION: (A-A)*U MOTIF.  PKA OF A7~6.1 AND IS               
REMARK 800 PARTIALLY PROTONATED AT PH 6.8, A7 NH1 MAY BE STABILIZED             
REMARK 800 BY U16 O2 (MODELS 7-12).                                             
REMARK 999                                                                      
REMARK 999 SEQUENCE                                                             
REMARK 999 CORRESPONDS TO GB D10027, GROUP II RNA BASES 1733   1753             
DBREF  17RA      1    21  PDB    17RA     17RA             1     21             
SEQRES   1     21    G   G   C   G   U   A   A   G   G   A   U   U   A          
SEQRES   2     21    C   C   U   A   U   G   C   C                              
SITE     1 AAU  3   A     6    A     7    U    16                               
CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1          
ORIGX1      1.000000  0.000000  0.000000        0.00000                         
ORIGX2      0.000000  1.000000  0.000000        0.00000                         
ORIGX3      0.000000  0.000000  1.000000        0.00000                         
SCALE1      1.000000  0.000000  0.000000        0.00000                         
SCALE2      0.000000  1.000000  0.000000        0.00000                         
SCALE3      0.000000  0.000000  1.000000        0.00000                         
MODEL        1                                                                  
ATOM      1  O5*   G     1      35.729  -4.043  -8.815  1.00  0.00           O  
ATOM      2  C5*   G     1      36.486  -5.197  -8.419  1.00  0.00           C  
ATOM      3  C4*   G     1      35.947  -6.491  -9.034  1.00  0.00           C  
ATOM      4  O4*   G     1      36.433  -7.650  -8.339  1.00  0.00           O  
ATOM      5  C3*   G     1      34.424  -6.588  -8.918  1.00  0.00           C  
ATOM      6  O3*   G     1      33.771  -5.934 -10.022  1.00  0.00           O  
ATOM      7  C2*   G     1      34.203  -8.095  -8.942  1.00  0.00           C  
ATOM      8  O2*   G     1      34.122  -8.584 -10.288  1.00  0.00           O  
ATOM      9  C1*   G     1      35.423  -8.665  -8.218  1.00  0.00           C  
ATOM     10  N9    G     1      35.126  -8.972  -6.792  1.00  0.00           N  
ATOM     11  C8    G     1      35.282  -8.197  -5.684  1.00  0.00           C  
ATOM     12  N7    G     1      34.953  -8.707  -4.545  1.00  0.00           N  
ATOM     13  C5    G     1      34.521  -9.985  -4.914  1.00  0.00           C  
ATOM     14  C6    G     1      34.026 -11.046  -4.105  1.00  0.00           C  
ATOM     15  O6    G     1      33.873 -11.066  -2.885  1.00  0.00           O  
ATOM     16  N1    G     1      33.701 -12.162  -4.869  1.00  0.00           N  
ATOM     17  C2    G     1      33.833 -12.251  -6.243  1.00  0.00           C  
ATOM     18  N2    G     1      33.467 -13.405  -6.801  1.00  0.00           N  
ATOM     19  N3    G     1      34.297 -11.256  -7.009  1.00  0.00           N  
ATOM     20  C4    G     1      34.621 -10.158  -6.287  1.00  0.00           C  
ATOM     21 1H5*   G     1      36.434  -5.308  -7.345  1.00  0.00           H  
ATOM     22 2H5*   G     1      37.533  -5.043  -8.723  1.00  0.00           H  
ATOM     23  H4*   G     1      36.243  -6.549 -10.082  1.00  0.00           H  
ATOM     24  H3*   G     1      34.089  -6.173  -7.960  1.00  0.00           H  
ATOM     25 1H2*   G     1      33.294  -8.345  -8.391  1.00  0.00           H  
ATOM     26 2HO*   G     1      33.856  -9.506 -10.242  1.00  0.00           H  
ATOM     27  H1*   G     1      35.757  -9.572  -8.727  1.00  0.00           H  
ATOM     28  H8    G     1      35.670  -7.179  -5.754  1.00  0.00           H  
ATOM     29  H1    G     1      33.342 -12.958  -4.361  1.00  0.00           H  
ATOM     30 1H2    G     1      33.541 -13.530  -7.800  1.00  0.00           H  
ATOM     31 2H2    G     1      33.116 -14.156  -6.224  1.00  0.00           H  
ATOM     32  H5T   G     1      35.846  -3.936  -9.762  1.00  0.00           H  
ATOM     33  P     G     2      32.263  -5.363  -9.892  1.00  0.00           P  
ATOM     34  O1P   G     2      31.976  -4.551 -11.096  1.00  0.00           O  
ATOM     35  O2P   G     2      32.103  -4.764  -8.548  1.00  0.00           O  
ATOM     36  O5*   G     2      31.365  -6.705  -9.956  1.00  0.00           O  
ATOM     37  C5*   G     2      31.011  -7.307 -11.214  1.00  0.00           C  
ATOM     38  C4*   G     2      30.343  -8.676 -11.050  1.00  0.00           C  
ATOM     39  O4*   G     2      31.035  -9.497 -10.097  1.00  0.00           O  
ATOM     40  C3*   G     2      28.923  -8.561 -10.491  1.00  0.00           C  
ATOM     41  O3*   G     2      27.951  -8.310 -11.523  1.00  0.00           O  
ATOM     42  C2*   G     2      28.731  -9.940  -9.871  1.00  0.00           C  
ATOM     43  O2*   G     2      28.290 -10.893 -10.847  1.00  0.00           O  
ATOM     44  C1*   G     2      30.119 -10.301  -9.334  1.00  0.00           C  
ATOM     45  N9    G     2      30.233 -10.024  -7.876  1.00  0.00           N  
ATOM     46  C8    G     2      30.688  -8.915  -7.234  1.00  0.00           C  
ATOM     47  N7    G     2      30.687  -8.928  -5.942  1.00  0.00           N  
ATOM     48  C5    G     2      30.171 -10.196  -5.665  1.00  0.00           C  
ATOM     49  C6    G     2      29.919 -10.822  -4.414  1.00  0.00           C  
ATOM     50  O6    G     2      30.109 -10.370  -3.286  1.00  0.00           O  
ATOM     51  N1    G     2      29.397 -12.100  -4.577  1.00  0.00           N  
ATOM     52  C2    G     2      29.145 -12.706  -5.793  1.00  0.00           C  
ATOM     53  N2    G     2      28.642 -13.941  -5.748  1.00  0.00           N  
ATOM     54  N3    G     2      29.379 -12.124  -6.974  1.00  0.00           N  
ATOM     55  C4    G     2      29.889 -10.878  -6.841  1.00  0.00           C  
ATOM     56 1H5*   G     2      31.908  -7.466 -11.801  1.00  0.00           H  
ATOM     57 2H5*   G     2      30.337  -6.621 -11.752  1.00  0.00           H  
ATOM     58  H4*   G     2      30.320  -9.189 -12.015  1.00  0.00           H  
ATOM     59  H3*   G     2      28.882  -7.791  -9.713  1.00  0.00           H  
ATOM     60 1H2*   G     2      28.016  -9.880  -9.045  1.00  0.00           H  
ATOM     61 2HO*   G     2      27.375 -10.684 -11.057  1.00  0.00           H  
ATOM     62  H1*   G     2      30.322 -11.357  -9.525  1.00  0.00           H  
ATOM     63  H8    G     2      31.045  -8.046  -7.786  1.00  0.00           H  
ATOM     64  H1    G     2      29.191 -12.607  -3.725  1.00  0.00           H  
ATOM     65 1H2    G     2      28.437 -14.436  -6.605  1.00  0.00           H  
ATOM     66 2H2    G     2      28.464 -14.382  -4.857  1.00  0.00           H  
ATOM     67  P     C     3      26.771  -7.224 -11.322  1.00  0.00           P  
ATOM     68  O1P   C     3      26.262  -6.847 -12.661  1.00  0.00           O  
ATOM     69  O2P   C     3      27.254  -6.170 -10.400  1.00  0.00           O  
ATOM     70  O5*   C     3      25.628  -8.075 -10.559  1.00  0.00           O  
ATOM     71  C5*   C     3      24.477  -8.572 -11.257  1.00  0.00           C  
ATOM     72  C4*   C     3      23.641  -9.523 -10.386  1.00  0.00           C  
ATOM     73  O4*   C     3      24.471 -10.338  -9.542  1.00  0.00           O  
ATOM     74  C3*   C     3      22.733  -8.770  -9.408  1.00  0.00           C  
ATOM     75  O3*   C     3      21.495  -8.361 -10.017  1.00  0.00           O  
ATOM     76  C2*   C     3      22.506  -9.840  -8.347  1.00  0.00           C  
ATOM     77  O2*   C     3      21.435 -10.719  -8.716  1.00  0.00           O  
ATOM     78  C1*   C     3      23.839 -10.592  -8.275  1.00  0.00           C  
ATOM     79  N1    C     3      24.678 -10.136  -7.128  1.00  0.00           N  
ATOM     80  C2    C     3      24.508 -10.789  -5.908  1.00  0.00           C  
ATOM     81  O2    C     3      23.685 -11.698  -5.800  1.00  0.00           O  
ATOM     82  N3    C     3      25.268 -10.395  -4.851  1.00  0.00           N  
ATOM     83  C4    C     3      26.159  -9.403  -4.971  1.00  0.00           C  
ATOM     84  N4    C     3      26.883  -9.048  -3.910  1.00  0.00           N  
ATOM     85  C5    C     3      26.343  -8.725  -6.217  1.00  0.00           C  
ATOM     86  C6    C     3      25.588  -9.118  -7.263  1.00  0.00           C  
ATOM     87 1H5*   C     3      24.807  -9.110 -12.148  1.00  0.00           H  
ATOM     88 2H5*   C     3      23.856  -7.730 -11.565  1.00  0.00           H  
ATOM     89  H4*   C     3      23.038 -10.170 -11.027  1.00  0.00           H  
ATOM     90  H3*   C     3      23.265  -7.915  -8.975  1.00  0.00           H  
ATOM     91 1H2*   C     3      22.294  -9.367  -7.384  1.00  0.00           H  
ATOM     92 2HO*   C     3      20.613 -10.251  -8.558  1.00  0.00           H  
ATOM     93  H1*   C     3      23.644 -11.663  -8.184  1.00  0.00           H  
ATOM     94 1H4    C     3      27.559  -8.301  -3.985  1.00  0.00           H  
ATOM     95 2H4    C     3      26.757  -9.525  -3.028  1.00  0.00           H  
ATOM     96  H5    C     3      27.069  -7.917  -6.317  1.00  0.00           H  
ATOM     97  H6    C     3      25.706  -8.619  -8.224  1.00  0.00           H  
ATOM     98  P     G     4      20.708  -7.037  -9.523  1.00  0.00           P  
ATOM     99  O1P   G     4      19.779  -6.624 -10.600  1.00  0.00           O  
ATOM    100  O2P   G     4      21.698  -6.069  -8.998  1.00  0.00           O  
ATOM    101  O5*   G     4      19.832  -7.583  -8.278  1.00  0.00           O  
ATOM    102  C5*   G     4      18.633  -8.350  -8.491  1.00  0.00           C  
ATOM    103  C4*   G     4      18.228  -9.169  -7.262  1.00  0.00           C  
ATOM    104  O4*   G     4      19.372  -9.689  -6.570  1.00  0.00           O  
ATOM    105  C3*   G     4      17.500  -8.321  -6.213  1.00  0.00           C  
ATOM    106  O3*   G     4      16.092  -8.212  -6.500  1.00  0.00           O  
ATOM    107  C2*   G     4      17.750  -9.137  -4.951  1.00  0.00           C  
ATOM    108  O2*   G     4      16.781 -10.184  -4.808  1.00  0.00           O  
ATOM    109  C1*   G     4      19.155  -9.706  -5.149  1.00  0.00           C  
ATOM    110  N9    G     4      20.181  -8.913  -4.420  1.00  0.00           N  
ATOM    111  C8    G     4      21.010  -7.930  -4.866  1.00  0.00           C  
ATOM    112  N7    G     4      21.826  -7.403  -4.014  1.00  0.00           N  
ATOM    113  C5    G     4      21.518  -8.108  -2.845  1.00  0.00           C  
ATOM    114  C6    G     4      22.071  -7.993  -1.538  1.00  0.00           C  
ATOM    115  O6    G     4      22.961  -7.238  -1.149  1.00  0.00           O  
ATOM    116  N1    G     4      21.477  -8.887  -0.653  1.00  0.00           N  
ATOM    117  C2    G     4      20.476  -9.781  -0.978  1.00  0.00           C  
ATOM    118  N2    G     4      20.036 -10.563   0.009  1.00  0.00           N  
ATOM    119  N3    G     4      19.952  -9.894  -2.202  1.00  0.00           N  
ATOM    120  C4    G     4      20.513  -9.033  -3.084  1.00  0.00           C  
ATOM    121 1H5*   G     4      18.796  -9.055  -9.296  1.00  0.00           H  
ATOM    122 2H5*   G     4      17.821  -7.658  -8.769  1.00  0.00           H  
ATOM    123  H4*   G     4      17.590 -10.000  -7.571  1.00  0.00           H  
ATOM    124  H3*   G     4      17.968  -7.333  -6.126  1.00  0.00           H  
ATOM    125 1H2*   G     4      17.734  -8.483  -4.076  1.00  0.00           H  
ATOM    126 2HO*   G     4      16.922 -10.590  -3.948  1.00  0.00           H  
ATOM    127  H1*   G     4      19.178 -10.735  -4.796  1.00  0.00           H  
ATOM    128  H8    G     4      20.986  -7.595  -5.904  1.00  0.00           H  
ATOM    129  H1    G     4      21.814  -8.859   0.298  1.00  0.00           H  
ATOM    130 1H2    G     4      20.432 -10.478   0.935  1.00  0.00           H  
ATOM    131 2H2    G     4      19.309 -11.239  -0.169  1.00  0.00           H  
ATOM    132  P     U     5      15.233  -6.916  -6.060  1.00  0.00           P  
ATOM    133  O1P   U     5      13.850  -7.096  -6.556  1.00  0.00           O  
ATOM    134  O2P   U     5      15.992  -5.701  -6.434  1.00  0.00           O  
ATOM    135  O5*   U     5      15.217  -7.025  -4.446  1.00  0.00           O  
ATOM    136  C5*   U     5      13.987  -7.142  -3.713  1.00  0.00           C  
ATOM    137  C4*   U     5      14.174  -6.810  -2.224  1.00  0.00           C  
ATOM    138  O4*   U     5      15.442  -7.271  -1.730  1.00  0.00           O  
ATOM    139  C3*   U     5      14.196  -5.299  -1.970  1.00  0.00           C  
ATOM    140  O3*   U     5      12.866  -4.773  -1.798  1.00  0.00           O  
ATOM    141  C2*   U     5      15.006  -5.214  -0.682  1.00  0.00           C  
ATOM    142  O2*   U     5      14.170  -5.389   0.469  1.00  0.00           O  
ATOM    143  C1*   U     5      16.027  -6.347  -0.796  1.00  0.00           C  
ATOM    144  N1    U     5      17.357  -5.858  -1.260  1.00  0.00           N  
ATOM    145  C2    U     5      18.197  -5.285  -0.314  1.00  0.00           C  
ATOM    146  O2    U     5      17.870  -5.156   0.866  1.00  0.00           O  
ATOM    147  N3    U     5      19.433  -4.863  -0.773  1.00  0.00           N  
ATOM    148  C4    U     5      19.901  -4.963  -2.073  1.00  0.00           C  
ATOM    149  O4    U     5      21.023  -4.555  -2.368  1.00  0.00           O  
ATOM    150  C5    U     5      18.963  -5.572  -2.991  1.00  0.00           C  
ATOM    151  C6    U     5      17.746  -5.989  -2.569  1.00  0.00           C  
ATOM    152 1H5*   U     5      13.620  -8.166  -3.803  1.00  0.00           H  
ATOM    153 2H5*   U     5      13.248  -6.464  -4.143  1.00  0.00           H  
ATOM    154  H4*   U     5      13.372  -7.270  -1.643  1.00  0.00           H  
ATOM    155  H3*   U     5      14.723  -4.785  -2.783  1.00  0.00           H  
ATOM    156 1H2*   U     5      15.521  -4.251  -0.634  1.00  0.00           H  
ATOM    157 2HO*   U     5      14.701  -5.187   1.242  1.00  0.00           H  
ATOM    158  H1*   U     5      16.136  -6.836   0.175  1.00  0.00           H  
ATOM    159  H3    U     5      20.051  -4.441  -0.096  1.00  0.00           H  
ATOM    160  H5    U     5      19.246  -5.709  -4.036  1.00  0.00           H  
ATOM    161  H6    U     5      17.062  -6.432  -3.288  1.00  0.00           H  
ATOM    162  P     A     6      12.565  -3.186  -1.896  1.00  0.00           P  
ATOM    163  O1P   A     6      11.096  -2.998  -1.909  1.00  0.00           O  
ATOM    164  O2P   A     6      13.385  -2.620  -2.991  1.00  0.00           O  
ATOM    165  O5*   A     6      13.139  -2.622  -0.494  1.00  0.00           O  
ATOM    166  C5*   A     6      12.398  -2.764   0.728  1.00  0.00           C  
ATOM    167  C4*   A     6      13.218  -2.316   1.947  1.00  0.00           C  
ATOM    168  O4*   A     6      14.583  -2.747   1.849  1.00  0.00           O  
ATOM    169  C3*   A     6      13.286  -0.789   2.053  1.00  0.00           C  
ATOM    170  O3*   A     6      12.314  -0.331   3.011  1.00  0.00           O  
ATOM    171  C2*   A     6      14.715  -0.482   2.497  1.00  0.00           C  
ATOM    172  O2*   A     6      14.736   0.082   3.814  1.00  0.00           O  
ATOM    173  C1*   A     6      15.474  -1.806   2.465  1.00  0.00           C  
ATOM    174  N9    A     6      16.741  -1.676   1.704  1.00  0.00           N  
ATOM    175  C8    A     6      16.951  -1.695   0.363  1.00  0.00           C  
ATOM    176  N7    A     6      18.166  -1.554  -0.053  1.00  0.00           N  
ATOM    177  C5    A     6      18.866  -1.420   1.150  1.00  0.00           C  
ATOM    178  C6    A     6      20.222  -1.233   1.445  1.00  0.00           C  
ATOM    179  N6    A     6      21.168  -1.143   0.507  1.00  0.00           N  
ATOM    180  N1    A     6      20.569  -1.141   2.740  1.00  0.00           N  
ATOM    181  C2    A     6      19.641  -1.226   3.696  1.00  0.00           C  
ATOM    182  N3    A     6      18.333  -1.403   3.532  1.00  0.00           N  
ATOM    183  C4    A     6      18.008  -1.492   2.225  1.00  0.00           C  
ATOM    184 1H5*   A     6      12.119  -3.812   0.856  1.00  0.00           H  
ATOM    185 2H5*   A     6      11.489  -2.162   0.668  1.00  0.00           H  
ATOM    186  H4*   A     6      12.776  -2.729   2.856  1.00  0.00           H  
ATOM    187  H3*   A     6      13.099  -0.339   1.074  1.00  0.00           H  
ATOM    188 1H2*   A     6      15.172   0.211   1.787  1.00  0.00           H  
ATOM    189 2HO*   A     6      14.285  -0.534   4.397  1.00  0.00           H  
ATOM    190  H1*   A     6      15.694  -2.126   3.485  1.00  0.00           H  
ATOM    191  H8    A     6      16.123  -1.811  -0.335  1.00  0.00           H  
ATOM    192 1H6    A     6      22.132  -1.008   0.774  1.00  0.00           H  
ATOM    193 2H6    A     6      20.918  -1.210  -0.469  1.00  0.00           H  
ATOM    194  H2    A     6      19.993  -1.143   4.724  1.00  0.00           H  
ATOM    195  P     A     7      11.860   1.219   3.079  1.00  0.00           P  
ATOM    196  O1P   A     7      10.531   1.280   3.729  1.00  0.00           O  
ATOM    197  O2P   A     7      12.051   1.826   1.741  1.00  0.00           O  
ATOM    198  O5*   A     7      12.943   1.858   4.095  1.00  0.00           O  
ATOM    199  C5*   A     7      13.656   3.061   3.772  1.00  0.00           C  
ATOM    200  C4*   A     7      14.906   3.224   4.652  1.00  0.00           C  
ATOM    201  O4*   A     7      15.883   2.209   4.374  1.00  0.00           O  
ATOM    202  C3*   A     7      15.639   4.540   4.375  1.00  0.00           C  
ATOM    203  O3*   A     7      15.112   5.625   5.160  1.00  0.00           O  
ATOM    204  C2*   A     7      17.066   4.195   4.787  1.00  0.00           C  
ATOM    205  O2*   A     7      17.264   4.395   6.194  1.00  0.00           O  
ATOM    206  C1*   A     7      17.224   2.721   4.409  1.00  0.00           C  
ATOM    207  N9    A     7      17.905   2.560   3.096  1.00  0.00           N  
ATOM    208  C8    A     7      17.377   2.501   1.843  1.00  0.00           C  
ATOM    209  N7    A     7      18.197   2.371   0.854  1.00  0.00           N  
ATOM    210  C5    A     7      19.433   2.334   1.510  1.00  0.00           C  
ATOM    211  C6    A     7      20.748   2.208   1.046  1.00  0.00           C  
ATOM    212  N6    A     7      21.058   2.096  -0.246  1.00  0.00           N  
ATOM    213  N1    A     7      21.731   2.204   1.964  1.00  0.00           N  
ATOM    214  C2    A     7      21.445   2.319   3.262  1.00  0.00           C  
ATOM    215  N3    A     7      20.239   2.446   3.811  1.00  0.00           N  
ATOM    216  C4    A     7      19.268   2.446   2.873  1.00  0.00           C  
ATOM    217 1H5*   A     7      12.996   3.918   3.924  1.00  0.00           H  
ATOM    218 2H5*   A     7      13.959   3.027   2.725  1.00  0.00           H  
ATOM    219  H4*   A     7      14.622   3.171   5.705  1.00  0.00           H  
ATOM    220  H3*   A     7      15.604   4.774   3.305  1.00  0.00           H  
ATOM    221 1H2*   A     7      17.771   4.801   4.215  1.00  0.00           H  
ATOM    222 2HO*   A     7      16.670   3.797   6.655  1.00  0.00           H  
ATOM    223  H1*   A     7      17.791   2.203   5.185  1.00  0.00           H  
ATOM    224  H8    A     7      16.301   2.551   1.679  1.00  0.00           H  
ATOM    225 1H6    A     7      22.024   2.006  -0.527  1.00  0.00           H  
ATOM    226 2H6    A     7      20.327   2.100  -0.943  1.00  0.00           H  
ATOM    227  H2    A     7      22.293   2.303   3.948  1.00  0.00           H  
ATOM    228  P     G     8      15.175   7.152   4.632  1.00  0.00           P  
ATOM    229  O1P   G     8      14.224   7.956   5.434  1.00  0.00           O  
ATOM    230  O2P   G     8      15.064   7.145   3.156  1.00  0.00           O  
ATOM    231  O5*   G     8      16.679   7.602   5.021  1.00  0.00           O  
ATOM    232  C5*   G     8      16.974   8.190   6.301  1.00  0.00           C  
ATOM    233  C4*   G     8      18.459   8.097   6.664  1.00  0.00           C  
ATOM    234  O4*   G     8      19.038   6.861   6.220  1.00  0.00           O  
ATOM    235  C3*   G     8      19.294   9.178   5.970  1.00  0.00           C  
ATOM    236  O3*   G     8      19.301  10.432   6.677  1.00  0.00           O  
ATOM    237  C2*   G     8      20.676   8.536   5.995  1.00  0.00           C  
ATOM    238  O2*   G     8      21.341   8.783   7.242  1.00  0.00           O  
ATOM    239  C1*   G     8      20.400   7.045   5.800  1.00  0.00           C  
ATOM    240  N9    G     8      20.603   6.633   4.384  1.00  0.00           N  
ATOM    241  C8    G     8      19.716   6.588   3.353  1.00  0.00           C  
ATOM    242  N7    G     8      20.153   6.196   2.203  1.00  0.00           N  
ATOM    243  C5    G     8      21.500   5.939   2.479  1.00  0.00           C  
ATOM    244  C6    G     8      22.532   5.475   1.617  1.00  0.00           C  
ATOM    245  O6    G     8      22.457   5.193   0.423  1.00  0.00           O  
ATOM    246  N1    G     8      23.743   5.352   2.291  1.00  0.00           N  
ATOM    247  C2    G     8      23.940   5.638   3.629  1.00  0.00           C  
ATOM    248  N2    G     8      25.175   5.455   4.101  1.00  0.00           N  
ATOM    249  N3    G     8      22.975   6.075   4.444  1.00  0.00           N  
ATOM    250  C4    G     8      21.786   6.203   3.810  1.00  0.00           C  
ATOM    251 1H5*   G     8      16.430   7.661   7.074  1.00  0.00           H  
ATOM    252 2H5*   G     8      16.653   9.244   6.283  1.00  0.00           H  
ATOM    253  H4*   G     8      18.576   8.179   7.746  1.00  0.00           H  
ATOM    254  H3*   G     8      18.960   9.308   4.933  1.00  0.00           H  
ATOM    255 1H2*   G     8      21.274   8.913   5.164  1.00  0.00           H  
ATOM    256 2HO*   G     8      22.250   8.489   7.145  1.00  0.00           H  
ATOM    257  H1*   G     8      21.060   6.464   6.449  1.00  0.00           H  
ATOM    258  H8    G     8      18.671   6.869   3.491  1.00  0.00           H  
ATOM    259  H1    G     8      24.528   5.030   1.739  1.00  0.00           H  
ATOM    260 1H2    G     8      25.377   5.646   5.072  1.00  0.00           H  
ATOM    261 2H2    G     8      25.906   5.124   3.488  1.00  0.00           H  
ATOM    262  P     G     9      19.503  11.838   5.902  1.00  0.00           P  
ATOM    263  O1P   G     9      19.258  12.934   6.867  1.00  0.00           O  
ATOM    264  O2P   G     9      18.733  11.791   4.637  1.00  0.00           O  
ATOM    265  O5*   G     9      21.078  11.827   5.533  1.00  0.00           O  
ATOM    266  C5*   G     9      22.079  12.117   6.527  1.00  0.00           C  
ATOM    267  C4*   G     9      23.469  11.607   6.136  1.00  0.00           C  
ATOM    268  O4*   G     9      23.399  10.353   5.441  1.00  0.00           O  
ATOM    269  C3*   G     9      24.176  12.539   5.154  1.00  0.00           C  
ATOM    270  O3*   G     9      24.856  13.639   5.786  1.00  0.00           O  
ATOM    271  C2*   G     9      25.176  11.572   4.531  1.00  0.00           C  
ATOM    272  O2*   G     9      26.330  11.407   5.367  1.00  0.00           O  
ATOM    273  C1*   G     9      24.390  10.265   4.404  1.00  0.00           C  
ATOM    274  N9    G     9      23.764  10.121   3.058  1.00  0.00           N  
ATOM    275  C8    G     9      22.457  10.216   2.689  1.00  0.00           C  
ATOM    276  N7    G     9      22.173  10.051   1.440  1.00  0.00           N  
ATOM    277  C5    G     9      23.433   9.814   0.884  1.00  0.00           C  
ATOM    278  C6    G     9      23.801   9.558  -0.468  1.00  0.00           C  
ATOM    279  O6    G     9      23.073   9.498  -1.460  1.00  0.00           O  
ATOM    280  N1    G     9      25.171   9.376  -0.598  1.00  0.00           N  
ATOM    281  C2    G     9      26.084   9.433   0.440  1.00  0.00           C  
ATOM    282  N2    G     9      27.363   9.241   0.117  1.00  0.00           N  
ATOM    283  N3    G     9      25.747   9.673   1.710  1.00  0.00           N  
ATOM    284  C4    G     9      24.415   9.854   1.864  1.00  0.00           C  
ATOM    285 1H5*   G     9      21.816  11.626   7.453  1.00  0.00           H  
ATOM    286 2H5*   G     9      22.106  13.208   6.687  1.00  0.00           H  
ATOM    287  H4*   G     9      24.079  11.487   7.029  1.00  0.00           H  
ATOM    288  H3*   G     9      23.471  12.891   4.392  1.00  0.00           H  
ATOM    289 1H2*   G     9      25.472  11.929   3.543  1.00  0.00           H  
ATOM    290 2HO*   G     9      26.061  10.888   6.129  1.00  0.00           H  
ATOM    291  H1*   G     9      25.050   9.416   4.600  1.00  0.00           H  
ATOM    292  H8    G     9      21.675  10.422   3.420  1.00  0.00           H  
ATOM    293  H1    G     9      25.505   9.184  -1.532  1.00  0.00           H  
ATOM    294 1H2    G     9      28.075   9.277   0.832  1.00  0.00           H  
ATOM    295 2H2    G     9      27.620   9.060  -0.844  1.00  0.00           H  
ATOM    296  P     A    10      24.574  15.170   5.346  1.00  0.00           P  
ATOM    297  O1P   A    10      25.088  16.058   6.414  1.00  0.00           O  
ATOM    298  O2P   A    10      23.160  15.284   4.920  1.00  0.00           O  
ATOM    299  O5*   A    10      25.509  15.351   4.040  1.00  0.00           O  
ATOM    300  C5*   A    10      26.926  15.558   4.158  1.00  0.00           C  
ATOM    301  C4*   A    10      27.695  14.893   3.010  1.00  0.00           C  
ATOM    302  O4*   A    10      27.088  13.644   2.639  1.00  0.00           O  
ATOM    303  C3*   A    10      27.666  15.734   1.730  1.00  0.00           C  
ATOM    304  O3*   A    10      28.738  16.681   1.557  1.00  0.00           O  
ATOM    305  C2*   A    10      27.848  14.655   0.676  1.00  0.00           C  
ATOM    306  O2*   A    10      29.235  14.386   0.433  1.00  0.00           O  
ATOM    307  C1*   A    10      27.132  13.433   1.221  1.00  0.00           C  
ATOM    308  N9    A    10      25.782  13.325   0.619  1.00  0.00           N  
ATOM    309  C8    A    10      24.567  13.233   1.208  1.00  0.00           C  
ATOM    310  N7    A    10      23.531  13.159   0.438  1.00  0.00           N  
ATOM    311  C5    A    10      24.118  13.212  -0.830  1.00  0.00           C  
ATOM    312  C6    A    10      23.590  13.185  -2.128  1.00  0.00           C  
ATOM    313  N6    A    10      22.285  13.087  -2.383  1.00  0.00           N  
ATOM    314  N1    A    10      24.460  13.266  -3.151  1.00  0.00           N  
ATOM    315  C2    A    10      25.772  13.368  -2.917  1.00  0.00           C  
ATOM    316  N3    A    10      26.373  13.394  -1.738  1.00  0.00           N  
ATOM    317  C4    A    10      25.485  13.314  -0.728  1.00  0.00           C  
ATOM    318 1H5*   A    10      27.271  15.137   5.101  1.00  0.00           H  
ATOM    319 2H5*   A    10      27.130  16.630   4.154  1.00  0.00           H  
ATOM    320  H4*   A    10      28.728  14.715   3.315  1.00  0.00           H  
ATOM    321  H3*   A    10      26.687  16.215   1.611  1.00  0.00           H  
ATOM    322 1H2*   A    10      27.356  14.965  -0.242  1.00  0.00           H  
ATOM    323 2HO*   A    10      29.611  14.051   1.250  1.00  0.00           H  
ATOM    324  H1*   A    10      27.712  12.533   0.999  1.00  0.00           H  
ATOM    325  H8    A    10      24.467  13.253   2.288  1.00  0.00           H  
ATOM    326 1H6    A    10      21.957  13.072  -3.339  1.00  0.00           H  
ATOM    327 2H6    A    10      21.624  13.026  -1.622  1.00  0.00           H  
ATOM    328  H2    A    10      26.420  13.467  -3.792  1.00  0.00           H  
ATOM    329  P     U    11      28.576  17.995   0.621  1.00  0.00           P  
ATOM    330  O1P   U    11      29.925  18.403   0.166  1.00  0.00           O  
ATOM    331  O2P   U    11      27.725  18.967   1.343  1.00  0.00           O  
ATOM    332  O5*   U    11      27.742  17.456  -0.660  1.00  0.00           O  
ATOM    333  C5*   U    11      28.394  16.929  -1.831  1.00  0.00           C  
ATOM    334  C4*   U    11      27.827  17.520  -3.137  1.00  0.00           C  
ATOM    335  O4*   U    11      26.520  17.016  -3.444  1.00  0.00           O  
ATOM    336  C3*   U    11      27.619  19.039  -3.065  1.00  0.00           C  
ATOM    337  O3*   U    11      27.694  19.750  -4.306  1.00  0.00           O  
ATOM    338  C2*   U    11      26.178  19.201  -2.604  1.00  0.00           C  
ATOM    339  O2*   U    11      25.622  20.454  -3.029  1.00  0.00           O  
ATOM    340  C1*   U    11      25.499  18.016  -3.297  1.00  0.00           C  
ATOM    341  N1    U    11      24.332  17.511  -2.520  1.00  0.00           N  
ATOM    342  C2    U    11      23.089  17.546  -3.137  1.00  0.00           C  
ATOM    343  O2    U    11      22.939  17.962  -4.286  1.00  0.00           O  
ATOM    344  N3    U    11      22.021  17.082  -2.388  1.00  0.00           N  
ATOM    345  C4    U    11      22.082  16.594  -1.094  1.00  0.00           C  
ATOM    346  O4    U    11      21.063  16.209  -0.521  1.00  0.00           O  
ATOM    347  C5    U    11      23.411  16.592  -0.526  1.00  0.00           C  
ATOM    348  C6    U    11      24.475  17.039  -1.237  1.00  0.00           C  
ATOM    349 1H5*   U    11      28.260  15.848  -1.848  1.00  0.00           H  
ATOM    350 2H5*   U    11      29.459  17.145  -1.780  1.00  0.00           H  
ATOM    351  H4*   U    11      28.499  17.279  -3.968  1.00  0.00           H  
ATOM    352  H3*   U    11      28.293  19.482  -2.333  1.00  0.00           H  
ATOM    353 1H2*   U    11      26.115  19.093  -1.518  1.00  0.00           H  
ATOM    354 2HO*   U    11      25.589  20.443  -3.989  1.00  0.00           H  
ATOM    355  H1*   U    11      25.169  18.324  -4.292  1.00  0.00           H  
ATOM    356  H3    U    11      21.111  17.104  -2.827  1.00  0.00           H  
ATOM    357  H5    U    11      23.559  16.227   0.491  1.00  0.00           H  
ATOM    358  H6    U    11      25.464  17.021  -0.780  1.00  0.00           H  
ATOM    359  P     U    12      29.104  20.030  -5.029  1.00  0.00           P  
ATOM    360  O1P   U    12      30.182  19.371  -4.256  1.00  0.00           O  
ATOM    361  O2P   U    12      29.190  21.477  -5.328  1.00  0.00           O  
ATOM    362  O5*   U    12      28.904  19.242  -6.418  1.00  0.00           O  
ATOM    363  C5*   U    12      29.559  17.999  -6.687  1.00  0.00           C  
ATOM    364  C4*   U    12      28.656  17.042  -7.476  1.00  0.00           C  
ATOM    365  O4*   U    12      27.374  17.616  -7.774  1.00  0.00           O  
ATOM    366  C3*   U    12      29.275  16.676  -8.821  1.00  0.00           C  
ATOM    367  O3*   U    12      29.008  15.408  -9.437  1.00  0.00           O  
ATOM    368  C2*   U    12      28.653  17.722  -9.741  1.00  0.00           C  
ATOM    369  O2*   U    12      28.547  17.237 -11.086  1.00  0.00           O  
ATOM    370  C1*   U    12      27.274  18.004  -9.154  1.00  0.00           C  
ATOM    371  N1    U    12      26.895  19.438  -9.299  1.00  0.00           N  
ATOM    372  C2    U    12      26.259  19.820 -10.473  1.00  0.00           C  
ATOM    373  O2    U    12      26.014  19.022 -11.377  1.00  0.00           O  
ATOM    374  N3    U    12      25.916  21.158 -10.572  1.00  0.00           N  
ATOM    375  C4    U    12      26.148  22.135  -9.619  1.00  0.00           C  
ATOM    376  O4    U    12      25.799  23.298  -9.814  1.00  0.00           O  
ATOM    377  C5    U    12      26.815  21.652  -8.430  1.00  0.00           C  
ATOM    378  C6    U    12      27.161  20.347  -8.306  1.00  0.00           C  
ATOM    379 1H5*   U    12      29.833  17.528  -5.743  1.00  0.00           H  
ATOM    380 2H5*   U    12      30.462  18.191  -7.254  1.00  0.00           H  
ATOM    381  H4*   U    12      28.507  16.133  -6.892  1.00  0.00           H  
ATOM    382  H3*   U    12      30.355  16.846  -8.759  1.00  0.00           H  
ATOM    383 1H2*   U    12      29.243  18.624  -9.718  1.00  0.00           H  
ATOM    384 2HO*   U    12      28.308  17.983 -11.641  1.00  0.00           H  
ATOM    385  H1*   U    12      26.532  17.378  -9.652  1.00  0.00           H  
ATOM    386  H3    U    12      25.452  21.448 -11.422  1.00  0.00           H  
ATOM    387  H5    U    12      27.041  22.348  -7.621  1.00  0.00           H  
ATOM    388  H6    U    12      27.668  20.016  -7.400  1.00  0.00           H  
ATOM    389  P     A    13      29.931  14.883 -10.657  1.00  0.00           P  
ATOM    390  O1P   A    13      30.926  15.931 -10.983  1.00  0.00           O  
ATOM    391  O2P   A    13      29.042  14.363 -11.721  1.00  0.00           O  
ATOM    392  O5*   A    13      30.717  13.645  -9.995  1.00  0.00           O  
ATOM    393  C5*   A    13      30.316  12.299 -10.252  1.00  0.00           C  
ATOM    394  C4*   A    13      29.638  11.672  -9.027  1.00  0.00           C  
ATOM    395  O4*   A    13      28.904  12.629  -8.255  1.00  0.00           O  
ATOM    396  C3*   A    13      28.602  10.612  -9.409  1.00  0.00           C  
ATOM    397  O3*   A    13      28.539   9.603  -8.396  1.00  0.00           O  
ATOM    398  C2*   A    13      27.285  11.373  -9.478  1.00  0.00           C  
ATOM    399  O2*   A    13      26.163  10.519  -9.217  1.00  0.00           O  
ATOM    400  C1*   A    13      27.494  12.395  -8.358  1.00  0.00           C  
ATOM    401  N9    A    13      26.739  13.655  -8.572  1.00  0.00           N  
ATOM    402  C8    A    13      26.342  14.278  -9.713  1.00  0.00           C  
ATOM    403  N7    A    13      25.576  15.311  -9.592  1.00  0.00           N  
ATOM    404  C5    A    13      25.456  15.411  -8.203  1.00  0.00           C  
ATOM    405  C6    A    13      24.783  16.309  -7.365  1.00  0.00           C  
ATOM    406  N6    A    13      24.048  17.324  -7.820  1.00  0.00           N  
ATOM    407  N1    A    13      24.896  16.115  -6.041  1.00  0.00           N  
ATOM    408  C2    A    13      25.627  15.104  -5.571  1.00  0.00           C  
ATOM    409  N3    A    13      26.290  14.202  -6.269  1.00  0.00           N  
ATOM    410  C4    A    13      26.166  14.412  -7.584  1.00  0.00           C  
ATOM    411 1H5*   A    13      31.198  11.713 -10.511  1.00  0.00           H  
ATOM    412 2H5*   A    13      29.625  12.287 -11.095  1.00  0.00           H  
ATOM    413  H4*   A    13      30.393  11.224  -8.384  1.00  0.00           H  
ATOM    414  H3*   A    13      28.839  10.174 -10.382  1.00  0.00           H  
ATOM    415 1H2*   A    13      27.187  11.875 -10.443  1.00  0.00           H  
ATOM    416 2HO*   A    13      26.029   9.977  -9.998  1.00  0.00           H  
ATOM    417  H1*   A    13      27.154  11.948  -7.421  1.00  0.00           H  
ATOM    418  H8    A    13      26.726  13.980 -10.681  1.00  0.00           H  
ATOM    419 1H6    A    13      23.582  17.942  -7.171  1.00  0.00           H  
ATOM    420 2H6    A    13      23.957  17.476  -8.815  1.00  0.00           H  
ATOM    421  H2    A    13      25.718  15.025  -4.491  1.00  0.00           H  
ATOM    422  P     C    14      28.746   8.040  -8.738  1.00  0.00           P  
ATOM    423  O1P   C    14      29.625   7.923  -9.924  1.00  0.00           O  
ATOM    424  O2P   C    14      27.417   7.388  -8.737  1.00  0.00           O  
ATOM    425  O5*   C    14      29.561   7.526  -7.444  1.00  0.00           O  
ATOM    426  C5*   C    14      30.941   7.870  -7.256  1.00  0.00           C  
ATOM    427  C4*   C    14      31.185   8.628  -5.942  1.00  0.00           C  
ATOM    428  O4*   C    14      30.361   9.798  -5.843  1.00  0.00           O  
ATOM    429  C3*   C    14      30.808   7.790  -4.720  1.00  0.00           C  
ATOM    430  O3*   C    14      31.916   6.964  -4.319  1.00  0.00           O  
ATOM    431  C2*   C    14      30.486   8.862  -3.686  1.00  0.00           C  
ATOM    432  O2*   C    14      31.671   9.303  -3.009  1.00  0.00           O  
ATOM    433  C1*   C    14      29.864   9.996  -4.509  1.00  0.00           C  
ATOM    434  N1    C    14      28.377   9.960  -4.487  1.00  0.00           N  
ATOM    435  C2    C    14      27.730  10.276  -3.297  1.00  0.00           C  
ATOM    436  O2    C    14      28.380  10.549  -2.289  1.00  0.00           O  
ATOM    437  N3    C    14      26.370  10.282  -3.285  1.00  0.00           N  
ATOM    438  C4    C    14      25.667  10.000  -4.388  1.00  0.00           C  
ATOM    439  N4    C    14      24.334   9.999  -4.323  1.00  0.00           N  
ATOM    440  C5    C    14      26.326   9.681  -5.621  1.00  0.00           C  
ATOM    441  C6    C    14      27.674   9.666  -5.621  1.00  0.00           C  
ATOM    442 1H5*   C    14      31.274   8.490  -8.090  1.00  0.00           H  
ATOM    443 2H5*   C    14      31.522   6.955  -7.247  1.00  0.00           H  
ATOM    444  H4*   C    14      32.234   8.923  -5.879  1.00  0.00           H  
ATOM    445  H3*   C    14      29.919   7.188  -4.929  1.00  0.00           H  
ATOM    446 1H2*   C    14      29.758   8.477  -2.969  1.00  0.00           H  
ATOM    447 2HO*   C    14      31.394   9.875  -2.291  1.00  0.00           H  
ATOM    448  H1*   C    14      30.211  10.959  -4.133  1.00  0.00           H  
ATOM    449 1H4    C    14      23.788   9.774  -5.140  1.00  0.00           H  
ATOM    450 2H4    C    14      23.870  10.216  -3.449  1.00  0.00           H  
ATOM    451  H5    C    14      25.761   9.508  -6.542  1.00  0.00           H  
ATOM    452  H6    C    14      28.212   9.392  -6.529  1.00  0.00           H  
ATOM    453  P     C    15      31.719   5.402  -3.957  1.00  0.00           P  
ATOM    454  O1P   C    15      33.013   4.715  -4.167  1.00  0.00           O  
ATOM    455  O2P   C    15      30.505   4.908  -4.646  1.00  0.00           O  
ATOM    456  O5*   C    15      31.415   5.441  -2.371  1.00  0.00           O  
ATOM    457  C5*   C    15      32.446   5.757  -1.422  1.00  0.00           C  
ATOM    458  C4*   C    15      31.871   6.019  -0.022  1.00  0.00           C  
ATOM    459  O4*   C    15      30.766   6.936  -0.065  1.00  0.00           O  
ATOM    460  C3*   C    15      31.284   4.752   0.604  1.00  0.00           C  
ATOM    461  O3*   C    15      32.277   3.986   1.310  1.00  0.00           O  
ATOM    462  C2*   C    15      30.252   5.327   1.564  1.00  0.00           C  
ATOM    463  O2*   C    15      30.848   5.659   2.826  1.00  0.00           O  
ATOM    464  C1*   C    15      29.722   6.578   0.856  1.00  0.00           C  
ATOM    465  N1    C    15      28.433   6.316   0.152  1.00  0.00           N  
ATOM    466  C2    C    15      27.307   6.054   0.933  1.00  0.00           C  
ATOM    467  O2    C    15      27.389   6.058   2.160  1.00  0.00           O  
ATOM    468  N3    C    15      26.125   5.807   0.305  1.00  0.00           N  
ATOM    469  C4    C    15      26.039   5.814  -1.032  1.00  0.00           C  
ATOM    470  N4    C    15      24.864   5.567  -1.612  1.00  0.00           N  
ATOM    471  C5    C    15      27.188   6.085  -1.843  1.00  0.00           C  
ATOM    472  C6    C    15      28.356   6.330  -1.216  1.00  0.00           C  
ATOM    473 1H5*   C    15      32.980   6.647  -1.757  1.00  0.00           H  
ATOM    474 2H5*   C    15      33.149   4.923  -1.368  1.00  0.00           H  
ATOM    475  H4*   C    15      32.651   6.426   0.625  1.00  0.00           H  
ATOM    476  H3*   C    15      30.788   4.145  -0.162  1.00  0.00           H  
ATOM    477 1H2*   C    15      29.443   4.608   1.705  1.00  0.00           H  
ATOM    478 2HO*   C    15      30.135   5.881   3.429  1.00  0.00           H  
ATOM    479  H1*   C    15      29.590   7.384   1.583  1.00  0.00           H  
ATOM    480 1H4    C    15      24.787   5.570  -2.619  1.00  0.00           H  
ATOM    481 2H4    C    15      24.052   5.375  -1.044  1.00  0.00           H  
ATOM    482  H5    C    15      27.120   6.093  -2.932  1.00  0.00           H  
ATOM    483  H6    C    15      29.245   6.553  -1.807  1.00  0.00           H  
ATOM    484  P     U    16      32.239   2.370   1.327  1.00  0.00           P  
ATOM    485  O1P   U    16      33.340   1.894   2.196  1.00  0.00           O  
ATOM    486  O2P   U    16      32.151   1.891  -0.072  1.00  0.00           O  
ATOM    487  O5*   U    16      30.835   2.048   2.060  1.00  0.00           O  
ATOM    488  C5*   U    16      30.717   2.070   3.491  1.00  0.00           C  
ATOM    489  C4*   U    16      29.251   1.992   3.949  1.00  0.00           C  
ATOM    490  O4*   U    16      28.384   2.751   3.092  1.00  0.00           O  
ATOM    491  C3*   U    16      28.697   0.566   3.878  1.00  0.00           C  
ATOM    492  O3*   U    16      29.025  -0.175   5.067  1.00  0.00           O  
ATOM    493  C2*   U    16      27.200   0.829   3.756  1.00  0.00           C  
ATOM    494  O2*   U    16      26.599   1.011   5.045  1.00  0.00           O  
ATOM    495  C1*   U    16      27.110   2.108   2.919  1.00  0.00           C  
ATOM    496  N1    U    16      26.823   1.816   1.483  1.00  0.00           N  
ATOM    497  C2    U    16      25.498   1.595   1.127  1.00  0.00           C  
ATOM    498  O2    U    16      24.584   1.626   1.950  1.00  0.00           O  
ATOM    499  N3    U    16      25.262   1.334  -0.212  1.00  0.00           N  
ATOM    500  C4    U    16      26.213   1.273  -1.216  1.00  0.00           C  
ATOM    501  O4    U    16      25.887   1.032  -2.377  1.00  0.00           O  
ATOM    502  C5    U    16      27.564   1.513  -0.761  1.00  0.00           C  
ATOM    503  C6    U    16      27.827   1.773   0.545  1.00  0.00           C  
ATOM    504 1H5*   U    16      31.156   2.996   3.868  1.00  0.00           H  
ATOM    505 2H5*   U    16      31.268   1.225   3.908  1.00  0.00           H  
ATOM    506  H4*   U    16      29.167   2.371   4.969  1.00  0.00           H  
ATOM    507  H3*   U    16      29.068   0.058   2.979  1.00  0.00           H  
ATOM    508 1H2*   U    16      26.719   0.003   3.226  1.00  0.00           H  
ATOM    509 2HO*   U    16      25.648   1.034   4.917  1.00  0.00           H  
ATOM    510  H1*   U    16      26.326   2.752   3.326  1.00  0.00           H  
ATOM    511  H3    U    16      24.300   1.173  -0.483  1.00  0.00           H  
ATOM    512  H5    U    16      28.385   1.485  -1.479  1.00  0.00           H  
ATOM    513  H6    U    16      28.857   1.950   0.851  1.00  0.00           H  
ATOM    514  P     A    17      28.993  -1.792   5.089  1.00  0.00           P  
ATOM    515  O1P   A    17      29.952  -2.257   6.118  1.00  0.00           O  
ATOM    516  O2P   A    17      29.116  -2.284   3.697  1.00  0.00           O  
ATOM    517  O5*   A    17      27.495  -2.113   5.604  1.00  0.00           O  
ATOM    518  C5*   A    17      27.111  -1.892   6.970  1.00  0.00           C  
ATOM    519  C4*   A    17      25.653  -2.303   7.229  1.00  0.00           C  
ATOM    520  O4*   A    17      24.758  -1.701   6.281  1.00  0.00           O  
ATOM    521  C3*   A    17      25.445  -3.806   7.077  1.00  0.00           C  
ATOM    522  O3*   A    17      25.756  -4.478   8.311  1.00  0.00           O  
ATOM    523  C2*   A    17      23.962  -3.866   6.726  1.00  0.00           C  
ATOM    524  O2*   A    17      23.150  -3.848   7.907  1.00  0.00           O  
ATOM    525  C1*   A    17      23.728  -2.616   5.875  1.00  0.00           C  
ATOM    526  N9    A    17      23.815  -2.914   4.421  1.00  0.00           N  
ATOM    527  C8    A    17      24.831  -2.680   3.545  1.00  0.00           C  
ATOM    528  N7    A    17      24.655  -3.039   2.318  1.00  0.00           N  
ATOM    529  C5    A    17      23.368  -3.581   2.358  1.00  0.00           C  
ATOM    530  C6    A    17      22.551  -4.157   1.377  1.00  0.00           C  
ATOM    531  N6    A    17      22.927  -4.293   0.105  1.00  0.00           N  
ATOM    532  N1    A    17      21.335  -4.587   1.757  1.00  0.00           N  
ATOM    533  C2    A    17      20.941  -4.463   3.027  1.00  0.00           C  
ATOM    534  N3    A    17      21.631  -3.936   4.035  1.00  0.00           N  
ATOM    535  C4    A    17      22.846  -3.510   3.631  1.00  0.00           C  
ATOM    536 1H5*   A    17      27.226  -0.832   7.205  1.00  0.00           H  
ATOM    537 2H5*   A    17      27.766  -2.471   7.625  1.00  0.00           H  
ATOM    538  H4*   A    17      25.363  -2.005   8.231  1.00  0.00           H  
ATOM    539  H3*   A    17      26.048  -4.193   6.248  1.00  0.00           H  
ATOM    540 1H2*   A    17      23.755  -4.760   6.138  1.00  0.00           H  
ATOM    541 2HO*   A    17      22.244  -4.013   7.635  1.00  0.00           H  
ATOM    542  H1*   A    17      22.748  -2.193   6.109  1.00  0.00           H  
ATOM    543  H8    A    17      25.758  -2.205   3.865  1.00  0.00           H  
ATOM    544 1H6    A    17      22.300  -4.715  -0.565  1.00  0.00           H  
ATOM    545 2H6    A    17      23.839  -3.974  -0.191  1.00  0.00           H  
ATOM    546  H2    A    17      19.938  -4.827   3.262  1.00  0.00           H  
ATOM    547  P     U    18      26.538  -5.892   8.327  1.00  0.00           P  
ATOM    548  O1P   U    18      26.683  -6.322   9.737  1.00  0.00           O  
ATOM    549  O2P   U    18      27.740  -5.772   7.470  1.00  0.00           O  
ATOM    550  O5*   U    18      25.495  -6.886   7.598  1.00  0.00           O  
ATOM    551  C5*   U    18      24.294  -7.317   8.256  1.00  0.00           C  
ATOM    552  C4*   U    18      23.245  -7.824   7.257  1.00  0.00           C  
ATOM    553  O4*   U    18      23.032  -6.889   6.188  1.00  0.00           O  
ATOM    554  C3*   U    18      23.694  -9.107   6.552  1.00  0.00           C  
ATOM    555  O3*   U    18      23.374 -10.289   7.311  1.00  0.00           O  
ATOM    556  C2*   U    18      22.884  -9.044   5.263  1.00  0.00           C  
ATOM    557  O2*   U    18      21.575  -9.599   5.447  1.00  0.00           O  
ATOM    558  C1*   U    18      22.811  -7.551   4.933  1.00  0.00           C  
ATOM    559  N1    U    18      23.820  -7.152   3.908  1.00  0.00           N  
ATOM    560  C2    U    18      23.539  -7.458   2.581  1.00  0.00           C  
ATOM    561  O2    U    18      22.503  -8.034   2.247  1.00  0.00           O  
ATOM    562  N3    U    18      24.489  -7.077   1.651  1.00  0.00           N  
ATOM    563  C4    U    18      25.681  -6.425   1.918  1.00  0.00           C  
ATOM    564  O4    U    18      26.454  -6.134   1.007  1.00  0.00           O  
ATOM    565  C5    U    18      25.897  -6.142   3.319  1.00  0.00           C  
ATOM    566  C6    U    18      24.983  -6.503   4.256  1.00  0.00           C  
ATOM    567 1H5*   U    18      23.874  -6.478   8.815  1.00  0.00           H  
ATOM    568 2H5*   U    18      24.540  -8.119   8.955  1.00  0.00           H  
ATOM    569  H4*   U    18      22.300  -7.993   7.774  1.00  0.00           H  
ATOM    570  H3*   U    18      24.765  -9.062   6.325  1.00  0.00           H  
ATOM    571 1H2*   U    18      23.413  -9.576   4.470  1.00  0.00           H  
ATOM    572 2HO*   U    18      21.670 -10.555   5.471  1.00  0.00           H  
ATOM    573  H1*   U    18      21.807  -7.309   4.574  1.00  0.00           H  
ATOM    574  H3    U    18      24.291  -7.293   0.681  1.00  0.00           H  
ATOM    575  H5    U    18      26.810  -5.631   3.629  1.00  0.00           H  
ATOM    576  H6    U    18      25.173  -6.258   5.304  1.00  0.00           H  
ATOM    577  P     G    19      24.346 -11.580   7.312  1.00  0.00           P  
ATOM    578  O1P   G    19      24.010 -12.406   8.495  1.00  0.00           O  
ATOM    579  O2P   G    19      25.739 -11.120   7.112  1.00  0.00           O  
ATOM    580  O5*   G    19      23.883 -12.376   5.982  1.00  0.00           O  
ATOM    581  C5*   G    19      23.221 -13.652   6.066  1.00  0.00           C  
ATOM    582  C4*   G    19      22.809 -14.197   4.694  1.00  0.00           C  
ATOM    583  O4*   G    19      22.387 -13.150   3.808  1.00  0.00           O  
ATOM    584  C3*   G    19      23.980 -14.852   3.954  1.00  0.00           C  
ATOM    585  O3*   G    19      24.213 -16.220   4.336  1.00  0.00           O  
ATOM    586  C2*   G    19      23.494 -14.774   2.513  1.00  0.00           C  
ATOM    587  O2*   G    19      22.622 -15.866   2.197  1.00  0.00           O  
ATOM    588  C1*   G    19      22.756 -13.437   2.447  1.00  0.00           C  
ATOM    589  N9    G    19      23.620 -12.369   1.877  1.00  0.00           N  
ATOM    590  C8    G    19      24.415 -11.464   2.507  1.00  0.00           C  
ATOM    591  N7    G    19      25.082 -10.639   1.773  1.00  0.00           N  
ATOM    592  C5    G    19      24.701 -11.024   0.485  1.00  0.00           C  
ATOM    593  C6    G    19      25.095 -10.496  -0.775  1.00  0.00           C  
ATOM    594  O6    G    19      25.871  -9.570  -1.004  1.00  0.00           O  
ATOM    595  N1    G    19      24.481 -11.170  -1.824  1.00  0.00           N  
ATOM    596  C2    G    19      23.594 -12.221  -1.685  1.00  0.00           C  
ATOM    597  N2    G    19      23.107 -12.739  -2.815  1.00  0.00           N  
ATOM    598  N3    G    19      23.218 -12.722  -0.504  1.00  0.00           N  
ATOM    599  C4    G    19      23.805 -12.081   0.536  1.00  0.00           C  
ATOM    600 1H5*   G    19      22.311 -13.546   6.644  1.00  0.00           H  
ATOM    601 2H5*   G    19      23.894 -14.363   6.569  1.00  0.00           H  
ATOM    602  H4*   G    19      21.997 -14.916   4.812  1.00  0.00           H  
ATOM    603  H3*   G    19      24.887 -14.248   4.076  1.00  0.00           H  
ATOM    604 1H2*   G    19      24.350 -14.760   1.833  1.00  0.00           H  
ATOM    605 2HO*   G    19      22.486 -15.860   1.246  1.00  0.00           H  
ATOM    606  H1*   G    19      21.854 -13.544   1.842  1.00  0.00           H  
ATOM    607  H8    G    19      24.490 -11.435   3.595  1.00  0.00           H  
ATOM    608  H1    G    19      24.716 -10.852  -2.751  1.00  0.00           H  
ATOM    609 1H2    G    19      22.456 -13.510  -2.777  1.00  0.00           H  
ATOM    610 2H2    G    19      23.392 -12.361  -3.708  1.00  0.00           H  
ATOM    611  P     C    20      25.686 -16.880   4.233  1.00  0.00           P  
ATOM    612  O1P   C    20      25.634 -18.216   4.870  1.00  0.00           O  
ATOM    613  O2P   C    20      26.678 -15.885   4.698  1.00  0.00           O  
ATOM    614  O5*   C    20      25.882 -17.081   2.640  1.00  0.00           O  
ATOM    615  C5*   C    20      25.335 -18.223   1.964  1.00  0.00           C  
ATOM    616  C4*   C    20      25.611 -18.194   0.451  1.00  0.00           C  
ATOM    617  O4*   C    20      25.236 -16.940  -0.144  1.00  0.00           O  
ATOM    618  C3*   C    20      27.103 -18.328   0.127  1.00  0.00           C  
ATOM    619  O3*   C    20      27.549 -19.698   0.164  1.00  0.00           O  
ATOM    620  C2*   C    20      27.121 -17.766  -1.290  1.00  0.00           C  
ATOM    621  O2*   C    20      26.757 -18.765  -2.252  1.00  0.00           O  
ATOM    622  C1*   C    20      26.097 -16.627  -1.254  1.00  0.00           C  
ATOM    623  N1    C    20      26.747 -15.292  -1.105  1.00  0.00           N  
ATOM    624  C2    C    20      27.208 -14.670  -2.264  1.00  0.00           C  
ATOM    625  O2    C    20      27.079 -15.220  -3.356  1.00  0.00           O  
ATOM    626  N3    C    20      27.800 -13.449  -2.153  1.00  0.00           N  
ATOM    627  C4    C    20      27.940 -12.856  -0.961  1.00  0.00           C  
ATOM    628  N4    C    20      28.526 -11.659  -0.895  1.00  0.00           N  
ATOM    629  C5    C    20      27.471 -13.484   0.235  1.00  0.00           C  
ATOM    630  C6    C    20      26.885 -14.693   0.120  1.00  0.00           C  
ATOM    631 1H5*   C    20      24.256 -18.244   2.124  1.00  0.00           H  
ATOM    632 2H5*   C    20      25.771 -19.130   2.387  1.00  0.00           H  
ATOM    633  H4*   C    20      25.056 -18.999  -0.034  1.00  0.00           H  
ATOM    634  H3*   C    20      27.698 -17.695   0.796  1.00  0.00           H  
ATOM    635 1H2*   C    20      28.114 -17.367  -1.514  1.00  0.00           H  
ATOM    636 2HO*   C    20      25.840 -18.999  -2.091  1.00  0.00           H  
ATOM    637  H1*   C    20      25.508 -16.645  -2.175  1.00  0.00           H  
ATOM    638 1H4    C    20      28.640 -11.200  -0.004  1.00  0.00           H  
ATOM    639 2H4    C    20      28.858 -11.212  -1.738  1.00  0.00           H  
ATOM    640  H5    C    20      27.584 -13.002   1.207  1.00  0.00           H  
ATOM    641  H6    C    20      26.515 -15.194   1.013  1.00  0.00           H  
ATOM    642  P     C    21      29.106 -20.075   0.387  1.00  0.00           P  
ATOM    643  O1P   C    21      29.192 -21.528   0.664  1.00  0.00           O  
ATOM    644  O2P   C    21      29.688 -19.111   1.347  1.00  0.00           O  
ATOM    645  O5*   C    21      29.757 -19.794  -1.065  1.00  0.00           O  
ATOM    646  C5*   C    21      29.562 -20.704  -2.159  1.00  0.00           C  
ATOM    647  C4*   C    21      29.991 -20.090  -3.501  1.00  0.00           C  
ATOM    648  O4*   C    21      29.684 -18.689  -3.574  1.00  0.00           O  
ATOM    649  C3*   C    21      31.505 -20.153  -3.707  1.00  0.00           C  
ATOM    650  O3*   C    21      31.915 -21.419  -4.241  1.00  0.00           O  
ATOM    651  C2*   C    21      31.720 -19.028  -4.710  1.00  0.00           C  
ATOM    652  O2*   C    21      31.499 -19.481  -6.052  1.00  0.00           O  
ATOM    653  C1*   C    21      30.691 -17.968  -4.305  1.00  0.00           C  
ATOM    654  N1    C    21      31.301 -16.871  -3.494  1.00  0.00           N  
ATOM    655  C2    C    21      31.856 -15.796  -4.186  1.00  0.00           C  
ATOM    656  O2    C    21      31.842 -15.776  -5.415  1.00  0.00           O  
ATOM    657  N3    C    21      32.408 -14.779  -3.469  1.00  0.00           N  
ATOM    658  C4    C    21      32.423 -14.807  -2.131  1.00  0.00           C  
ATOM    659  N4    C    21      32.976 -13.791  -1.465  1.00  0.00           N  
ATOM    660  C5    C    21      31.858 -15.906  -1.409  1.00  0.00           C  
ATOM    661  C6    C    21      31.311 -16.910  -2.123  1.00  0.00           C  
ATOM    662 1H5*   C    21      28.504 -20.970  -2.214  1.00  0.00           H  
ATOM    663 2H5*   C    21      30.143 -21.610  -1.978  1.00  0.00           H  
ATOM    664  H4*   C    21      29.485 -20.610  -4.318  1.00  0.00           H  
ATOM    665  H3*   C    21      32.032 -19.932  -2.770  1.00  0.00           H  
ATOM    666 1H2*   C    21      32.733 -18.626  -4.605  1.00  0.00           H  
ATOM    667 2HO*   C    21      31.753 -18.767  -6.642  1.00  0.00           H  
ATOM    668  H1*   C    21      30.237 -17.545  -5.205  1.00  0.00           H  
ATOM    669 1H4    C    21      32.994 -13.798  -0.455  1.00  0.00           H  
ATOM    670 2H4    C    21      33.376 -13.013  -1.971  1.00  0.00           H  
ATOM    671  H5    C    21      31.869 -15.929  -0.318  1.00  0.00           H  
ATOM    672  H6    C    21      30.873 -17.760  -1.601  1.00  0.00           H  
ATOM    673  H3T   C    21      32.875 -21.422  -4.267  1.00  0.00           H  
TER     674        C    21                                                      
ENDMDL                                                                          
MODEL        2                                                                  
ATOM    675  O5*   G     1      37.008 -12.325 -11.644  1.00  0.00           O  
ATOM    676  C5*   G     1      35.778 -11.791 -11.135  1.00  0.00           C  
ATOM    677  C4*   G     1      34.582 -12.698 -11.437  1.00  0.00           C  
ATOM    678  O4*   G     1      34.720 -13.979 -10.804  1.00  0.00           O  
ATOM    679  C3*   G     1      33.275 -12.133 -10.875  1.00  0.00           C  
ATOM    680  O3*   G     1      32.656 -11.221 -11.801  1.00  0.00           O  
ATOM    681  C2*   G     1      32.443 -13.395 -10.683  1.00  0.00           C  
ATOM    682  O2*   G     1      31.758 -13.753 -11.891  1.00  0.00           O  
ATOM    683  C1*   G     1      33.467 -14.461 -10.291  1.00  0.00           C  
ATOM    684  N9    G     1      33.522 -14.660  -8.818  1.00  0.00           N  
ATOM    685  C8    G     1      34.370 -14.124  -7.898  1.00  0.00           C  
ATOM    686  N7    G     1      34.209 -14.467  -6.663  1.00  0.00           N  
ATOM    687  C5    G     1      33.121 -15.342  -6.745  1.00  0.00           C  
ATOM    688  C6    G     1      32.453 -16.062  -5.716  1.00  0.00           C  
ATOM    689  O6    G     1      32.698 -16.070  -4.511  1.00  0.00           O  
ATOM    690  N1    G     1      31.408 -16.828  -6.223  1.00  0.00           N  
ATOM    691  C2    G     1      31.047 -16.895  -7.556  1.00  0.00           C  
ATOM    692  N2    G     1      30.014 -17.687  -7.852  1.00  0.00           N  
ATOM    693  N3    G     1      31.669 -16.220  -8.527  1.00  0.00           N  
ATOM    694  C4    G     1      32.692 -15.467  -8.059  1.00  0.00           C  
ATOM    695 1H5*   G     1      35.579 -10.837 -11.606  1.00  0.00           H  
ATOM    696 2H5*   G     1      35.887 -11.650 -10.048  1.00  0.00           H  
ATOM    697  H4*   G     1      34.492 -12.838 -12.516  1.00  0.00           H  
ATOM    698  H3*   G     1      33.455 -11.650  -9.907  1.00  0.00           H  
ATOM    699 1H2*   G     1      31.730 -13.248  -9.869  1.00  0.00           H  
ATOM    700 2HO*   G     1      30.920 -14.147 -11.638  1.00  0.00           H  
ATOM    701  H1*   G     1      33.214 -15.405 -10.780  1.00  0.00           H  
ATOM    702  H8    G     1      35.160 -13.430  -8.191  1.00  0.00           H  
ATOM    703  H1    G     1      30.887 -17.371  -5.548  1.00  0.00           H  
ATOM    704 1H2    G     1      29.705 -17.775  -8.809  1.00  0.00           H  
ATOM    705 2H2    G     1      29.543 -18.197  -7.120  1.00  0.00           H  
ATOM    706  H5T   G     1      36.972 -12.266 -12.601  1.00  0.00           H  
ATOM    707  P     G     2      31.597 -10.105 -11.305  1.00  0.00           P  
ATOM    708  O1P   G     2      31.385  -9.148 -12.416  1.00  0.00           O  
ATOM    709  O2P   G     2      32.027  -9.603  -9.980  1.00  0.00           O  
ATOM    710  O5*   G     2      30.241 -10.961 -11.103  1.00  0.00           O  
ATOM    711  C5*   G     2      29.403 -11.311 -12.220  1.00  0.00           C  
ATOM    712  C4*   G     2      28.254 -12.246 -11.831  1.00  0.00           C  
ATOM    713  O4*   G     2      28.692 -13.298 -10.958  1.00  0.00           O  
ATOM    714  C3*   G     2      27.167 -11.521 -11.033  1.00  0.00           C  
ATOM    715  O3*   G     2      26.224 -10.863 -11.900  1.00  0.00           O  
ATOM    716  C2*   G     2      26.527 -12.672 -10.268  1.00  0.00           C  
ATOM    717  O2*   G     2      25.523 -13.325 -11.057  1.00  0.00           O  
ATOM    718  C1*   G     2      27.697 -13.614  -9.968  1.00  0.00           C  
ATOM    719  N9    G     2      28.219 -13.421  -8.588  1.00  0.00           N  
ATOM    720  C8    G     2      29.283 -12.693  -8.154  1.00  0.00           C  
ATOM    721  N7    G     2      29.533 -12.693  -6.887  1.00  0.00           N  
ATOM    722  C5    G     2      28.519 -13.519  -6.392  1.00  0.00           C  
ATOM    723  C6    G     2      28.244 -13.918  -5.055  1.00  0.00           C  
ATOM    724  O6    G     2      28.853 -13.616  -4.030  1.00  0.00           O  
ATOM    725  N1    G     2      27.135 -14.753  -4.989  1.00  0.00           N  
ATOM    726  C2    G     2      26.378 -15.157  -6.074  1.00  0.00           C  
ATOM    727  N2    G     2      25.346 -15.961  -5.811  1.00  0.00           N  
ATOM    728  N3    G     2      26.631 -14.786  -7.334  1.00  0.00           N  
ATOM    729  C4    G     2      27.709 -13.971  -7.423  1.00  0.00           C  
ATOM    730 1H5*   G     2      29.996 -11.836 -12.961  1.00  0.00           H  
ATOM    731 2H5*   G     2      29.003 -10.384 -12.662  1.00  0.00           H  
ATOM    732  H4*   G     2      27.820 -12.686 -12.731  1.00  0.00           H  
ATOM    733  H3*   G     2      27.622 -10.811 -10.332  1.00  0.00           H  
ATOM    734 1H2*   G     2      26.097 -12.302  -9.333  1.00  0.00           H  
ATOM    735 2HO*   G     2      25.111 -13.992 -10.501  1.00  0.00           H  
ATOM    736  H1*   G     2      27.375 -14.649 -10.099  1.00  0.00           H  
ATOM    737  H8    G     2      29.908 -12.132  -8.849  1.00  0.00           H  
ATOM    738  H1    G     2      26.879 -15.079  -4.068  1.00  0.00           H  
ATOM    739 1H2    G     2      24.757 -16.290  -6.562  1.00  0.00           H  
ATOM    740 2H2    G     2      25.155 -16.243  -4.859  1.00  0.00           H  
ATOM    741  P     C     3      25.601  -9.419 -11.527  1.00  0.00           P  
ATOM    742  O1P   C     3      25.067  -8.814 -12.769  1.00  0.00           O  
ATOM    743  O2P   C     3      26.592  -8.675 -10.717  1.00  0.00           O  
ATOM    744  O5*   C     3      24.357  -9.804 -10.570  1.00  0.00           O  
ATOM    745  C5*   C     3      23.094 -10.207 -11.119  1.00  0.00           C  
ATOM    746  C4*   C     3      22.246 -10.979 -10.096  1.00  0.00           C  
ATOM    747  O4*   C     3      23.056 -11.826  -9.267  1.00  0.00           O  
ATOM    748  C3*   C     3      21.535 -10.046  -9.110  1.00  0.00           C  
ATOM    749  O3*   C     3      20.295  -9.540  -9.639  1.00  0.00           O  
ATOM    750  C2*   C     3      21.302 -10.992  -7.937  1.00  0.00           C  
ATOM    751  O2*   C     3      20.101 -11.756  -8.115  1.00  0.00           O  
ATOM    752  C1*   C     3      22.538 -11.897  -7.928  1.00  0.00           C  
ATOM    753  N1    C     3      23.542 -11.462  -6.911  1.00  0.00           N  
ATOM    754  C2    C     3      23.422 -11.988  -5.627  1.00  0.00           C  
ATOM    755  O2    C     3      22.515 -12.776  -5.364  1.00  0.00           O  
ATOM    756  N3    C     3      24.325 -11.608  -4.684  1.00  0.00           N  
ATOM    757  C4    C     3      25.310 -10.749  -4.974  1.00  0.00           C  
ATOM    758  N4    C     3      26.172 -10.399  -4.019  1.00  0.00           N  
ATOM    759  C5    C     3      25.444 -10.201  -6.290  1.00  0.00           C  
ATOM    760  C6    C     3      24.547 -10.580  -7.223  1.00  0.00           C  
ATOM    761 1H5*   C     3      23.270 -10.847 -11.985  1.00  0.00           H  
ATOM    762 2H5*   C     3      22.547  -9.319 -11.442  1.00  0.00           H  
ATOM    763  H4*   C     3      21.509 -11.590 -10.621  1.00  0.00           H  
ATOM    764  H3*   C     3      22.202  -9.229  -8.810  1.00  0.00           H  
ATOM    765 1H2*   C     3      21.254 -10.421  -7.007  1.00  0.00           H  
ATOM    766 2HO*   C     3      20.243 -12.346  -8.859  1.00  0.00           H  
ATOM    767  H1*   C     3      22.232 -12.925  -7.721  1.00  0.00           H  
ATOM    768 1H4    C     3      26.918  -9.750  -4.223  1.00  0.00           H  
ATOM    769 2H4    C     3      26.078 -10.784  -3.088  1.00  0.00           H  
ATOM    770  H5    C     3      26.246  -9.501  -6.530  1.00  0.00           H  
ATOM    771  H6    C     3      24.625 -10.180  -8.233  1.00  0.00           H  
ATOM    772  P     G     4      19.644  -8.163  -9.094  1.00  0.00           P  
ATOM    773  O1P   G     4      18.604  -7.733 -10.057  1.00  0.00           O  
ATOM    774  O2P   G     4      20.741  -7.233  -8.737  1.00  0.00           O  
ATOM    775  O5*   G     4      18.912  -8.625  -7.728  1.00  0.00           O  
ATOM    776  C5*   G     4      17.694  -9.392  -7.759  1.00  0.00           C  
ATOM    777  C4*   G     4      17.404 -10.094  -6.429  1.00  0.00           C  
ATOM    778  O4*   G     4      18.608 -10.564  -5.802  1.00  0.00           O  
ATOM    779  C3*   G     4      16.784  -9.148  -5.396  1.00  0.00           C  
ATOM    780  O3*   G     4      15.354  -9.046  -5.548  1.00  0.00           O  
ATOM    781  C2*   G     4      17.148  -9.851  -4.094  1.00  0.00           C  
ATOM    782  O2*   G     4      16.192 -10.871  -3.771  1.00  0.00           O  
ATOM    783  C1*   G     4      18.527 -10.452  -4.372  1.00  0.00           C  
ATOM    784  N9    G     4      19.622  -9.612  -3.818  1.00  0.00           N  
ATOM    785  C8    G     4      20.416  -8.692  -4.429  1.00  0.00           C  
ATOM    786  N7    G     4      21.319  -8.109  -3.711  1.00  0.00           N  
ATOM    787  C5    G     4      21.118  -8.699  -2.458  1.00  0.00           C  
ATOM    788  C6    G     4      21.800  -8.482  -1.229  1.00  0.00           C  
ATOM    789  O6    G     4      22.733  -7.713  -0.997  1.00  0.00           O  
ATOM    790  N1    G     4      21.285  -9.279  -0.212  1.00  0.00           N  
ATOM    791  C2    G     4      20.244 -10.175  -0.356  1.00  0.00           C  
ATOM    792  N2    G     4      19.892 -10.852   0.738  1.00  0.00           N  
ATOM    793  N3    G     4      19.600 -10.385  -1.506  1.00  0.00           N  
ATOM    794  C4    G     4      20.082  -9.620  -2.513  1.00  0.00           C  
ATOM    795 1H5*   G     4      17.777 -10.167  -8.510  1.00  0.00           H  
ATOM    796 2H5*   G     4      16.864  -8.717  -8.020  1.00  0.00           H  
ATOM    797  H4*   G     4      16.737 -10.941  -6.601  1.00  0.00           H  
ATOM    798  H3*   G     4      17.262  -8.163  -5.442  1.00  0.00           H  
ATOM    799 1H2*   G     4      17.216  -9.122  -3.284  1.00  0.00           H  
ATOM    800 2HO*   G     4      16.390 -11.175  -2.881  1.00  0.00           H  
ATOM    801  H1*   G     4      18.578 -11.447  -3.932  1.00  0.00           H  
ATOM    802  H8    G     4      20.294  -8.449  -5.485  1.00  0.00           H  
ATOM    803  H1    G     4      21.718  -9.178   0.693  1.00  0.00           H  
ATOM    804 1H2    G     4      20.379 -10.694   1.610  1.00  0.00           H  
ATOM    805 2H2    G     4      19.138 -11.524   0.697  1.00  0.00           H  
ATOM    806  P     U     5      14.561  -7.681  -5.200  1.00  0.00           P  
ATOM    807  O1P   U     5      13.126  -7.896  -5.497  1.00  0.00           O  
ATOM    808  O2P   U     5      15.277  -6.549  -5.832  1.00  0.00           O  
ATOM    809  O5*   U     5      14.740  -7.557  -3.598  1.00  0.00           O  
ATOM    810  C5*   U     5      13.691  -7.942  -2.694  1.00  0.00           C  
ATOM    811  C4*   U     5      13.954  -7.447  -1.263  1.00  0.00           C  
ATOM    812  O4*   U     5      15.270  -7.801  -0.808  1.00  0.00           O  
ATOM    813  C3*   U     5      13.921  -5.918  -1.168  1.00  0.00           C  
ATOM    814  O3*   U     5      12.578  -5.438  -0.963  1.00  0.00           O  
ATOM    815  C2*   U     5      14.806  -5.663   0.046  1.00  0.00           C  
ATOM    816  O2*   U     5      14.052  -5.744   1.264  1.00  0.00           O  
ATOM    817  C1*   U     5      15.869  -6.760  -0.017  1.00  0.00           C  
ATOM    818  N1    U     5      17.145  -6.276  -0.615  1.00  0.00           N  
ATOM    819  C2    U     5      18.017  -5.575   0.208  1.00  0.00           C  
ATOM    820  O2    U     5      17.752  -5.322   1.383  1.00  0.00           O  
ATOM    821  N3    U     5      19.208  -5.171  -0.369  1.00  0.00           N  
ATOM    822  C4    U     5      19.604  -5.402  -1.677  1.00  0.00           C  
ATOM    823  O4    U     5      20.693  -4.997  -2.083  1.00  0.00           O  
ATOM    824  C5    U     5      18.638  -6.137  -2.463  1.00  0.00           C  
ATOM    825  C6    U     5      17.463  -6.542  -1.923  1.00  0.00           C  
ATOM    826 1H5*   U     5      13.616  -9.031  -2.684  1.00  0.00           H  
ATOM    827 2H5*   U     5      12.744  -7.527  -3.047  1.00  0.00           H  
ATOM    828  H4*   U     5      13.212  -7.875  -0.585  1.00  0.00           H  
ATOM    829  H3*   U     5      14.368  -5.472  -2.063  1.00  0.00           H  
ATOM    830 1H2*   U     5      15.277  -4.682  -0.042  1.00  0.00           H  
ATOM    831 2HO*   U     5      13.710  -6.639   1.331  1.00  0.00           H  
ATOM    832  H1*   U     5      16.058  -7.137   0.990  1.00  0.00           H  
ATOM    833  H3    U     5      19.849  -4.657   0.218  1.00  0.00           H  
ATOM    834  H5    U     5      18.863  -6.384  -3.502  1.00  0.00           H  
ATOM    835  H6    U     5      16.755  -7.083  -2.546  1.00  0.00           H  
ATOM    836  P     A     6      12.187  -3.888  -1.211  1.00  0.00           P  
ATOM    837  O1P   A     6      10.715  -3.766  -1.112  1.00  0.00           O  
ATOM    838  O2P   A     6      12.880  -3.421  -2.433  1.00  0.00           O  
ATOM    839  O5*   A     6      12.848  -3.140   0.060  1.00  0.00           O  
ATOM    840  C5*   A     6      12.211  -3.149   1.348  1.00  0.00           C  
ATOM    841  C4*   A     6      13.092  -2.495   2.425  1.00  0.00           C  
ATOM    842  O4*   A     6      14.452  -2.942   2.339  1.00  0.00           O  
ATOM    843  C3*   A     6      13.164  -0.974   2.252  1.00  0.00           C  
ATOM    844  O3*   A     6      12.183  -0.340   3.092  1.00  0.00           O  
ATOM    845  C2*   A     6      14.587  -0.603   2.659  1.00  0.00           C  
ATOM    846  O2*   A     6      14.605   0.044   3.937  1.00  0.00           O  
ATOM    847  C1*   A     6      15.381  -1.908   2.702  1.00  0.00           C  
ATOM    848  N9    A     6      16.530  -1.858   1.766  1.00  0.00           N  
ATOM    849  C8    A     6      16.552  -2.015   0.421  1.00  0.00           C  
ATOM    850  N7    A     6      17.692  -1.914  -0.176  1.00  0.00           N  
ATOM    851  C5    A     6      18.551  -1.655   0.899  1.00  0.00           C  
ATOM    852  C6    A     6      19.931  -1.434   0.977  1.00  0.00           C  
ATOM    853  N6    A     6      20.734  -1.440  -0.089  1.00  0.00           N  
ATOM    854  N1    A     6      20.454  -1.205   2.195  1.00  0.00           N  
ATOM    855  C2    A     6      19.670  -1.195   3.275  1.00  0.00           C  
ATOM    856  N3    A     6      18.355  -1.391   3.316  1.00  0.00           N  
ATOM    857  C4    A     6      17.854  -1.618   2.084  1.00  0.00           C  
ATOM    858 1H5*   A     6      12.007  -4.181   1.637  1.00  0.00           H  
ATOM    859 2H5*   A     6      11.265  -2.607   1.283  1.00  0.00           H  
ATOM    860  H4*   A     6      12.697  -2.734   3.414  1.00  0.00           H  
ATOM    861  H3*   A     6      12.998  -0.710   1.203  1.00  0.00           H  
ATOM    862 1H2*   A     6      15.020   0.055   1.901  1.00  0.00           H  
ATOM    863 2HO*   A     6      14.170  -0.543   4.561  1.00  0.00           H  
ATOM    864  H1*   A     6      15.744  -2.083   3.716  1.00  0.00           H  
ATOM    865  H8    A     6      15.636  -2.206  -0.137  1.00  0.00           H  
ATOM    866 1H6    A     6      21.723  -1.274   0.023  1.00  0.00           H  
ATOM    867 2H6    A     6      20.351  -1.610  -1.008  1.00  0.00           H  
ATOM    868  H2    A     6      20.162  -1.003   4.230  1.00  0.00           H  
ATOM    869  P     A     7      11.764   1.207   2.880  1.00  0.00           P  
ATOM    870  O1P   A     7      10.422   1.404   3.472  1.00  0.00           O  
ATOM    871  O2P   A     7      12.003   1.572   1.465  1.00  0.00           O  
ATOM    872  O5*   A     7      12.838   1.986   3.804  1.00  0.00           O  
ATOM    873  C5*   A     7      13.666   3.034   3.277  1.00  0.00           C  
ATOM    874  C4*   A     7      14.900   3.262   4.164  1.00  0.00           C  
ATOM    875  O4*   A     7      15.877   2.226   3.983  1.00  0.00           O  
ATOM    876  C3*   A     7      15.645   4.552   3.803  1.00  0.00           C  
ATOM    877  O3*   A     7      15.097   5.691   4.495  1.00  0.00           O  
ATOM    878  C2*   A     7      17.058   4.234   4.276  1.00  0.00           C  
ATOM    879  O2*   A     7      17.221   4.541   5.668  1.00  0.00           O  
ATOM    880  C1*   A     7      17.219   2.735   4.019  1.00  0.00           C  
ATOM    881  N9    A     7      17.943   2.467   2.747  1.00  0.00           N  
ATOM    882  C8    A     7      17.459   2.277   1.491  1.00  0.00           C  
ATOM    883  N7    A     7      18.315   2.068   0.547  1.00  0.00           N  
ATOM    884  C5    A     7      19.527   2.120   1.243  1.00  0.00           C  
ATOM    885  C6    A     7      20.861   1.974   0.838  1.00  0.00           C  
ATOM    886  N6    A     7      21.219   1.739  -0.425  1.00  0.00           N  
ATOM    887  N1    A     7      21.811   2.082   1.784  1.00  0.00           N  
ATOM    888  C2    A     7      21.477   2.318   3.055  1.00  0.00           C  
ATOM    889  N3    A     7      20.250   2.473   3.547  1.00  0.00           N  
ATOM    890  C4    A     7      19.313   2.361   2.580  1.00  0.00           C  
ATOM    891 1H5*   A     7      13.083   3.955   3.226  1.00  0.00           H  
ATOM    892 2H5*   A     7      13.992   2.762   2.272  1.00  0.00           H  
ATOM    893  H4*   A     7      14.597   3.291   5.212  1.00  0.00           H  
ATOM    894  H3*   A     7      15.635   4.706   2.717  1.00  0.00           H  
ATOM    895 1H2*   A     7      17.782   4.792   3.676  1.00  0.00           H  
ATOM    896 2HO*   A     7      16.609   3.987   6.157  1.00  0.00           H  
ATOM    897  H1*   A     7      17.756   2.278   4.851  1.00  0.00           H  
ATOM    898  H8    A     7      16.387   2.290   1.287  1.00  0.00           H  
ATOM    899 1H6    A     7      22.196   1.643  -0.663  1.00  0.00           H  
ATOM    900 2H6    A     7      20.513   1.657  -1.143  1.00  0.00           H  
ATOM    901  H2    A     7      22.299   2.386   3.768  1.00  0.00           H  
ATOM    902  P     G     8      15.447   7.203   4.041  1.00  0.00           P  
ATOM    903  O1P   G     8      14.403   8.100   4.587  1.00  0.00           O  
ATOM    904  O2P   G     8      15.726   7.206   2.587  1.00  0.00           O  
ATOM    905  O5*   G     8      16.832   7.501   4.820  1.00  0.00           O  
ATOM    906  C5*   G     8      16.850   7.784   6.232  1.00  0.00           C  
ATOM    907  C4*   G     8      18.238   7.596   6.853  1.00  0.00           C  
ATOM    908  O4*   G     8      18.934   6.484   6.270  1.00  0.00           O  
ATOM    909  C3*   G     8      19.154   8.797   6.596  1.00  0.00           C  
ATOM    910  O3*   G     8      18.995   9.830   7.589  1.00  0.00           O  
ATOM    911  C2*   G     8      20.534   8.156   6.693  1.00  0.00           C  
ATOM    912  O2*   G     8      20.992   8.107   8.051  1.00  0.00           O  
ATOM    913  C1*   G     8      20.338   6.751   6.124  1.00  0.00           C  
ATOM    914  N9    G     8      20.771   6.671   4.703  1.00  0.00           N  
ATOM    915  C8    G     8      20.045   6.813   3.561  1.00  0.00           C  
ATOM    916  N7    G     8      20.673   6.707   2.436  1.00  0.00           N  
ATOM    917  C5    G     8      21.984   6.462   2.858  1.00  0.00           C  
ATOM    918  C6    G     8      23.163   6.253   2.090  1.00  0.00           C  
ATOM    919  O6    G     8      23.284   6.243   0.867  1.00  0.00           O  
ATOM    920  N1    G     8      24.271   6.042   2.904  1.00  0.00           N  
ATOM    921  C2    G     8      24.251   6.032   4.287  1.00  0.00           C  
ATOM    922  N2    G     8      25.415   5.807   4.898  1.00  0.00           N  
ATOM    923  N3    G     8      23.147   6.230   5.013  1.00  0.00           N  
ATOM    924  C4    G     8      22.054   6.438   4.243  1.00  0.00           C  
ATOM    925 1H5*   G     8      16.184   7.101   6.745  1.00  0.00           H  
ATOM    926 2H5*   G     8      16.503   8.818   6.385  1.00  0.00           H  
ATOM    927  H4*   G     8      18.138   7.432   7.928  1.00  0.00           H  
ATOM    928  H3*   G     8      18.991   9.190   5.586  1.00  0.00           H  
ATOM    929 1H2*   G     8      21.242   8.710   6.074  1.00  0.00           H  
ATOM    930 2HO*   G     8      20.386   7.545   8.538  1.00  0.00           H  
ATOM    931  H1*   G     8      20.907   6.036   6.721  1.00  0.00           H  
ATOM    932  H8    G     8      18.973   7.005   3.590  1.00  0.00           H  
ATOM    933  H1    G     8      25.150   5.891   2.427  1.00  0.00           H  
ATOM    934 1H2    G     8      25.460   5.785   5.906  1.00  0.00           H  
ATOM    935 2H2    G     8      26.250   5.655   4.351  1.00  0.00           H  
ATOM    936  P     G     9      19.371  11.371   7.273  1.00  0.00           P  
ATOM    937  O1P   G     9      18.904  12.198   8.410  1.00  0.00           O  
ATOM    938  O2P   G     9      18.919  11.692   5.901  1.00  0.00           O  
ATOM    939  O5*   G     9      20.988  11.363   7.282  1.00  0.00           O  
ATOM    940  C5*   G     9      21.729  11.337   8.516  1.00  0.00           C  
ATOM    941  C4*   G     9      23.200  10.954   8.316  1.00  0.00           C  
ATOM    942  O4*   G     9      23.358   9.961   7.292  1.00  0.00           O  
ATOM    943  C3*   G     9      24.044  12.132   7.832  1.00  0.00           C  
ATOM    944  O3*   G     9      24.504  12.988   8.895  1.00  0.00           O  
ATOM    945  C2*   G     9      25.211  11.387   7.193  1.00  0.00           C  
ATOM    946  O2*   G     9      26.160  10.963   8.181  1.00  0.00           O  
ATOM    947  C1*   G     9      24.545  10.187   6.512  1.00  0.00           C  
ATOM    948  N9    G     9      24.217  10.466   5.085  1.00  0.00           N  
ATOM    949  C8    G     9      23.010  10.681   4.496  1.00  0.00           C  
ATOM    950  N7    G     9      22.992  10.916   3.226  1.00  0.00           N  
ATOM    951  C5    G     9      24.352  10.856   2.902  1.00  0.00           C  
ATOM    952  C6    G     9      24.996  11.033   1.646  1.00  0.00           C  
ATOM    953  O6    G     9      24.487  11.290   0.555  1.00  0.00           O  
ATOM    954  N1    G     9      26.374  10.890   1.755  1.00  0.00           N  
ATOM    955  C2    G     9      27.054  10.612   2.927  1.00  0.00           C  
ATOM    956  N2    G     9      28.381  10.509   2.832  1.00  0.00           N  
ATOM    957  N3    G     9      26.456  10.446   4.111  1.00  0.00           N  
ATOM    958  C4    G     9      25.111  10.580   4.031  1.00  0.00           C  
ATOM    959 1H5*   G     9      21.300  10.595   9.177  1.00  0.00           H  
ATOM    960 2H5*   G     9      21.657  12.330   8.987  1.00  0.00           H  
ATOM    961  H4*   G     9      23.608  10.573   9.251  1.00  0.00           H  
ATOM    962  H3*   G     9      23.499  12.701   7.070  1.00  0.00           H  
ATOM    963 1H2*   G     9      25.694  12.023   6.448  1.00  0.00           H  
ATOM    964 2HO*   G     9      26.943  10.661   7.715  1.00  0.00           H  
ATOM    965  H1*   G     9      25.200   9.313   6.581  1.00  0.00           H  
ATOM    966  H8    G     9      22.084  10.649   5.073  1.00  0.00           H  
ATOM    967  H1    G     9      26.900  11.006   0.898  1.00  0.00           H  
ATOM    968 1H2    G     9      28.931  10.305   3.654  1.00  0.00           H  
ATOM    969 2H2    G     9      28.835  10.632   1.939  1.00  0.00           H  
ATOM    970  P     A    10      24.210  14.578   8.885  1.00  0.00           P  
ATOM    971  O1P   A    10      24.448  15.101  10.250  1.00  0.00           O  
ATOM    972  O2P   A    10      22.901  14.810   8.231  1.00  0.00           O  
ATOM    973  O5*   A    10      25.366  15.154   7.911  1.00  0.00           O  
ATOM    974  C5*   A    10      26.711  15.341   8.379  1.00  0.00           C  
ATOM    975  C4*   A    10      27.739  15.068   7.275  1.00  0.00           C  
ATOM    976  O4*   A    10      27.319  13.980   6.439  1.00  0.00           O  
ATOM    977  C3*   A    10      27.893  16.256   6.319  1.00  0.00           C  
ATOM    978  O3*   A    10      28.901  17.233   6.647  1.00  0.00           O  
ATOM    979  C2*   A    10      28.363  15.549   5.059  1.00  0.00           C  
ATOM    980  O2*   A    10      29.788  15.390   5.046  1.00  0.00           O  
ATOM    981  C1*   A    10      27.655  14.208   5.062  1.00  0.00           C  
ATOM    982  N9    A    10      26.458  14.265   4.188  1.00  0.00           N  
ATOM    983  C8    A    10      25.159  14.007   4.467  1.00  0.00           C  
ATOM    984  N7    A    10      24.302  14.150   3.511  1.00  0.00           N  
ATOM    985  C5    A    10      25.123  14.559   2.455  1.00  0.00           C  
ATOM    986  C6    A    10      24.862  14.897   1.121  1.00  0.00           C  
ATOM    987  N6    A    10      23.641  14.870   0.585  1.00  0.00           N  
ATOM    988  N1    A    10      25.910  15.265   0.361  1.00  0.00           N  
ATOM    989  C2    A    10      27.141  15.300   0.879  1.00  0.00           C  
ATOM    990  N3    A    10      27.497  14.996   2.117  1.00  0.00           N  
ATOM    991  C4    A    10      26.434  14.632   2.860  1.00  0.00           C  
ATOM    992 1H5*   A    10      26.899  14.662   9.209  1.00  0.00           H  
ATOM    993 2H5*   A    10      26.825  16.367   8.730  1.00  0.00           H  
ATOM    994  H4*   A    10      28.704  14.825   7.725  1.00  0.00           H  
ATOM    995  H3*   A    10      26.920  16.733   6.143  1.00  0.00           H  
ATOM    996 1H2*   A    10      28.038  16.114   4.190  1.00  0.00           H  
ATOM    997 2HO*   A    10      30.022  14.832   5.791  1.00  0.00           H  
ATOM    998  H1*   A    10      28.334  13.427   4.714  1.00  0.00           H  
ATOM    999  H8    A    10      24.848  13.720   5.467  1.00  0.00           H  
ATOM   1000 1H6    A    10      23.509  15.125  -0.384  1.00  0.00           H  
ATOM   1001 2H6    A    10      22.849  14.595   1.147  1.00  0.00           H  
ATOM   1002  H2    A    10      27.943  15.645   0.219  1.00  0.00           H  
ATOM   1003  P     U    11      28.812  18.764   6.118  1.00  0.00           P  
ATOM   1004  O1P   U    11      30.180  19.331   6.126  1.00  0.00           O  
ATOM   1005  O2P   U    11      27.736  19.447   6.871  1.00  0.00           O  
ATOM   1006  O5*   U    11      28.328  18.615   4.578  1.00  0.00           O  
ATOM   1007  C5*   U    11      29.263  18.459   3.494  1.00  0.00           C  
ATOM   1008  C4*   U    11      28.940  19.378   2.299  1.00  0.00           C  
ATOM   1009  O4*   U    11      27.773  18.954   1.580  1.00  0.00           O  
ATOM   1010  C3*   U    11      28.599  20.813   2.722  1.00  0.00           C  
ATOM   1011  O3*   U    11      28.865  21.847   1.767  1.00  0.00           O  
ATOM   1012  C2*   U    11      27.084  20.814   2.869  1.00  0.00           C  
ATOM   1013  O2*   U    11      26.535  22.124   2.678  1.00  0.00           O  
ATOM   1014  C1*   U    11      26.667  19.857   1.751  1.00  0.00           C  
ATOM   1015  N1    U    11      25.405  19.136   2.076  1.00  0.00           N  
ATOM   1016  C2    U    11      24.329  19.311   1.216  1.00  0.00           C  
ATOM   1017  O2    U    11      24.398  20.021   0.213  1.00  0.00           O  
ATOM   1018  N3    U    11      23.165  18.639   1.550  1.00  0.00           N  
ATOM   1019  C4    U    11      22.980  17.820   2.651  1.00  0.00           C  
ATOM   1020  O4    U    11      21.896  17.271   2.850  1.00  0.00           O  
ATOM   1021  C5    U    11      24.148  17.691   3.492  1.00  0.00           C  
ATOM   1022  C6    U    11      25.302  18.336   3.190  1.00  0.00           C  
ATOM   1023 1H5*   U    11      29.236  17.423   3.156  1.00  0.00           H  
ATOM   1024 2H5*   U    11      30.267  18.682   3.849  1.00  0.00           H  
ATOM   1025  H4*   U    11      29.789  19.390   1.610  1.00  0.00           H  
ATOM   1026  H3*   U    11      29.063  21.050   3.677  1.00  0.00           H  
ATOM   1027 1H2*   U    11      26.799  20.406   3.842  1.00  0.00           H  
ATOM   1028 2HO*   U    11      25.592  22.066   2.848  1.00  0.00           H  
ATOM   1029  H1*   U    11      26.534  20.425   0.826  1.00  0.00           H  
ATOM   1030  H3    U    11      22.374  18.759   0.931  1.00  0.00           H  
ATOM   1031  H5    U    11      24.099  17.067   4.386  1.00  0.00           H  
ATOM   1032  H6    U    11      26.165  18.216   3.844  1.00  0.00           H  
ATOM   1033  P     U    12      30.365  22.352   1.474  1.00  0.00           P  
ATOM   1034  O1P   U    12      31.310  21.554   2.289  1.00  0.00           O  
ATOM   1035  O2P   U    12      30.381  23.829   1.564  1.00  0.00           O  
ATOM   1036  O5*   U    12      30.538  21.941  -0.073  1.00  0.00           O  
ATOM   1037  C5*   U    12      31.344  20.830  -0.476  1.00  0.00           C  
ATOM   1038  C4*   U    12      30.748  20.111  -1.694  1.00  0.00           C  
ATOM   1039  O4*   U    12      29.536  20.722  -2.157  1.00  0.00           O  
ATOM   1040  C3*   U    12      31.705  20.135  -2.878  1.00  0.00           C  
ATOM   1041  O3*   U    12      31.712  19.032  -3.792  1.00  0.00           O  
ATOM   1042  C2*   U    12      31.243  21.382  -3.630  1.00  0.00           C  
ATOM   1043  O2*   U    12      31.513  21.273  -5.033  1.00  0.00           O  
ATOM   1044  C1*   U    12      29.743  21.469  -3.368  1.00  0.00           C  
ATOM   1045  N1    U    12      29.291  22.881  -3.225  1.00  0.00           N  
ATOM   1046  C2    U    12      28.934  23.562  -4.381  1.00  0.00           C  
ATOM   1047  O2    U    12      28.991  23.041  -5.495  1.00  0.00           O  
ATOM   1048  N3    U    12      28.512  24.870  -4.213  1.00  0.00           N  
ATOM   1049  C4    U    12      28.414  25.550  -3.010  1.00  0.00           C  
ATOM   1050  O4    U    12      28.024  26.716  -2.978  1.00  0.00           O  
ATOM   1051  C5    U    12      28.805  24.768  -1.857  1.00  0.00           C  
ATOM   1052  C6    U    12      29.223  23.486  -1.993  1.00  0.00           C  
ATOM   1053 1H5*   U    12      31.417  20.123   0.350  1.00  0.00           H  
ATOM   1054 2H5*   U    12      32.337  21.183  -0.719  1.00  0.00           H  
ATOM   1055  H4*   U    12      30.540  19.074  -1.425  1.00  0.00           H  
ATOM   1056  H3*   U    12      32.714  20.307  -2.500  1.00  0.00           H  
ATOM   1057 1H2*   U    12      31.729  22.258  -3.222  1.00  0.00           H  
ATOM   1058 2HO*   U    12      32.366  21.685  -5.193  1.00  0.00           H  
ATOM   1059  H1*   U    12      29.201  20.993  -4.184  1.00  0.00           H  
ATOM   1060  H3    U    12      28.249  25.377  -5.047  1.00  0.00           H  
ATOM   1061  H5    U    12      28.763  25.218  -0.864  1.00  0.00           H  
ATOM   1062  H6    U    12      29.517  22.924  -1.105  1.00  0.00           H  
ATOM   1063  P     A    13      33.089  18.602  -4.523  1.00  0.00           P  
ATOM   1064  O1P   A    13      33.745  17.561  -3.700  1.00  0.00           O  
ATOM   1065  O2P   A    13      33.835  19.830  -4.881  1.00  0.00           O  
ATOM   1066  O5*   A    13      32.559  17.931  -5.890  1.00  0.00           O  
ATOM   1067  C5*   A    13      32.848  16.569  -6.220  1.00  0.00           C  
ATOM   1068  C4*   A    13      32.210  15.588  -5.224  1.00  0.00           C  
ATOM   1069  O4*   A    13      31.454  16.258  -4.215  1.00  0.00           O  
ATOM   1070  C3*   A    13      31.229  14.633  -5.904  1.00  0.00           C  
ATOM   1071  O3*   A    13      31.260  13.370  -5.229  1.00  0.00           O  
ATOM   1072  C2*   A    13      29.888  15.340  -5.756  1.00  0.00           C  
ATOM   1073  O2*   A    13      28.797  14.411  -5.753  1.00  0.00           O  
ATOM   1074  C1*   A    13      30.062  16.013  -4.393  1.00  0.00           C  
ATOM   1075  N9    A    13      29.280  17.261  -4.260  1.00  0.00           N  
ATOM   1076  C8    A    13      29.030  18.262  -5.142  1.00  0.00           C  
ATOM   1077  N7    A    13      28.209  19.188  -4.770  1.00  0.00           N  
ATOM   1078  C5    A    13      27.885  18.771  -3.476  1.00  0.00           C  
ATOM   1079  C6    A    13      27.062  19.304  -2.476  1.00  0.00           C  
ATOM   1080  N6    A    13      26.356  20.425  -2.629  1.00  0.00           N  
ATOM   1081  N1    A    13      26.999  18.633  -1.315  1.00  0.00           N  
ATOM   1082  C2    A    13      27.702  17.509  -1.151  1.00  0.00           C  
ATOM   1083  N3    A    13      28.490  16.925  -2.024  1.00  0.00           N  
ATOM   1084  C4    A    13      28.541  17.608  -3.169  1.00  0.00           C  
ATOM   1085 1H5*   A    13      33.931  16.424  -6.218  1.00  0.00           H  
ATOM   1086 2H5*   A    13      32.471  16.364  -7.221  1.00  0.00           H  
ATOM   1087  H4*   A    13      32.986  15.009  -4.730  1.00  0.00           H  
ATOM   1088  H3*   A    13      31.481  14.517  -6.962  1.00  0.00           H  
ATOM   1089 1H2*   A    13      29.767  16.088  -6.540  1.00  0.00           H  
ATOM   1090 2HO*   A    13      27.985  14.922  -5.797  1.00  0.00           H  
ATOM   1091  H1*   A    13      29.736  15.316  -3.619  1.00  0.00           H  
ATOM   1092  H8    A    13      29.582  18.344  -6.068  1.00  0.00           H  
ATOM   1093 1H6    A    13      25.778  20.764  -1.872  1.00  0.00           H  
ATOM   1094 2H6    A    13      26.397  20.936  -3.499  1.00  0.00           H  
ATOM   1095  H2    A    13      27.656  17.039  -0.178  1.00  0.00           H  
ATOM   1096  P     C    14      31.030  11.984  -6.021  1.00  0.00           P  
ATOM   1097  O1P   C    14      31.972  11.929  -7.163  1.00  0.00           O  
ATOM   1098  O2P   C    14      29.579  11.813  -6.253  1.00  0.00           O  
ATOM   1099  O5*   C    14      31.496  10.901  -4.923  1.00  0.00           O  
ATOM   1100  C5*   C    14      32.862  10.835  -4.492  1.00  0.00           C  
ATOM   1101  C4*   C    14      33.022  11.139  -2.994  1.00  0.00           C  
ATOM   1102  O4*   C    14      32.307  12.320  -2.612  1.00  0.00           O  
ATOM   1103  C3*   C    14      32.430  10.031  -2.122  1.00  0.00           C  
ATOM   1104  O3*   C    14      33.398   8.990  -1.904  1.00  0.00           O  
ATOM   1105  C2*   C    14      32.093  10.783  -0.841  1.00  0.00           C  
ATOM   1106  O2*   C    14      33.228  10.854   0.032  1.00  0.00           O  
ATOM   1107  C1*   C    14      31.678  12.176  -1.328  1.00  0.00           C  
ATOM   1108  N1    C    14      30.202  12.319  -1.441  1.00  0.00           N  
ATOM   1109  C2    C    14      29.447  12.280  -0.272  1.00  0.00           C  
ATOM   1110  O2    C    14      29.995  12.119   0.819  1.00  0.00           O  
ATOM   1111  N3    C    14      28.101  12.439  -0.369  1.00  0.00           N  
ATOM   1112  C4    C    14      27.510  12.635  -1.555  1.00  0.00           C  
ATOM   1113  N4    C    14      26.183  12.760  -1.603  1.00  0.00           N  
ATOM   1114  C5    C    14      28.282  12.684  -2.763  1.00  0.00           C  
ATOM   1115  C6    C    14      29.615  12.514  -2.658  1.00  0.00           C  
ATOM   1116 1H5*   C    14      33.452  11.553  -5.065  1.00  0.00           H  
ATOM   1117 2H5*   C    14      33.235   9.837  -4.691  1.00  0.00           H  
ATOM   1118  H4*   C    14      34.081  11.268  -2.759  1.00  0.00           H  
ATOM   1119  H3*   C    14      31.517   9.634  -2.579  1.00  0.00           H  
ATOM   1120 1H2*   C    14      31.254  10.295  -0.339  1.00  0.00           H  
ATOM   1121 2HO*   C    14      32.922  11.211   0.870  1.00  0.00           H  
ATOM   1122  H1*   C    14      32.066  12.937  -0.649  1.00  0.00           H  
ATOM   1123 1H4    C    14      25.720  12.890  -2.489  1.00  0.00           H  
ATOM   1124 2H4    C    14      25.640  12.718  -0.751  1.00  0.00           H  
ATOM   1125  H5    C    14      27.815  12.905  -3.727  1.00  0.00           H  
ATOM   1126  H6    C    14      30.233  12.508  -3.558  1.00  0.00           H  
ATOM   1127  P     C    15      32.984   7.429  -1.947  1.00  0.00           P  
ATOM   1128  O1P   C    15      34.217   6.631  -2.132  1.00  0.00           O  
ATOM   1129  O2P   C    15      31.858   7.267  -2.896  1.00  0.00           O  
ATOM   1130  O5*   C    15      32.426   7.167  -0.454  1.00  0.00           O  
ATOM   1131  C5*   C    15      33.306   7.165   0.681  1.00  0.00           C  
ATOM   1132  C4*   C    15      32.536   7.006   2.001  1.00  0.00           C  
ATOM   1133  O4*   C    15      31.494   7.986   2.129  1.00  0.00           O  
ATOM   1134  C3*   C    15      31.811   5.661   2.084  1.00  0.00           C  
ATOM   1135  O3*   C    15      32.663   4.622   2.598  1.00  0.00           O  
ATOM   1136  C2*   C    15      30.681   5.984   3.053  1.00  0.00           C  
ATOM   1137  O2*   C    15      31.107   5.831   4.414  1.00  0.00           O  
ATOM   1138  C1*   C    15      30.317   7.442   2.750  1.00  0.00           C  
ATOM   1139  N1    C    15      29.122   7.547   1.862  1.00  0.00           N  
ATOM   1140  C2    C    15      27.889   7.166   2.390  1.00  0.00           C  
ATOM   1141  O2    C    15      27.809   6.774   3.552  1.00  0.00           O  
ATOM   1142  N3    C    15      26.789   7.250   1.594  1.00  0.00           N  
ATOM   1143  C4    C    15      26.881   7.689   0.333  1.00  0.00           C  
ATOM   1144  N4    C    15      25.779   7.756  -0.416  1.00  0.00           N  
ATOM   1145  C5    C    15      28.141   8.085  -0.219  1.00  0.00           C  
ATOM   1146  C6    C    15      29.230   8.000   0.573  1.00  0.00           C  
ATOM   1147 1H5*   C    15      33.858   8.105   0.706  1.00  0.00           H  
ATOM   1148 2H5*   C    15      34.015   6.341   0.581  1.00  0.00           H  
ATOM   1149  H4*   C    15      33.227   7.106   2.841  1.00  0.00           H  
ATOM   1150  H3*   C    15      31.400   5.389   1.104  1.00  0.00           H  
ATOM   1151 1H2*   C    15      29.826   5.334   2.848  1.00  0.00           H  
ATOM   1152 2HO*   C    15      30.330   5.925   4.969  1.00  0.00           H  
ATOM   1153  H1*   C    15      30.127   7.972   3.687  1.00  0.00           H  
ATOM   1154 1H4    C    15      25.835   8.084  -1.370  1.00  0.00           H  
ATOM   1155 2H4    C    15      24.888   7.476  -0.031  1.00  0.00           H  
ATOM   1156  H5    C    15      28.218   8.445  -1.247  1.00  0.00           H  
ATOM   1157  H6    C    15      30.201   8.306   0.184  1.00  0.00           H  
ATOM   1158  P     U    16      32.481   3.079   2.148  1.00  0.00           P  
ATOM   1159  O1P   U    16      33.543   2.281   2.801  1.00  0.00           O  
ATOM   1160  O2P   U    16      32.335   3.037   0.675  1.00  0.00           O  
ATOM   1161  O5*   U    16      31.061   2.672   2.807  1.00  0.00           O  
ATOM   1162  C5*   U    16      30.943   2.388   4.209  1.00  0.00           C  
ATOM   1163  C4*   U    16      29.476   2.232   4.643  1.00  0.00           C  
ATOM   1164  O4*   U    16      28.611   3.142   3.944  1.00  0.00           O  
ATOM   1165  C3*   U    16      28.916   0.848   4.304  1.00  0.00           C  
ATOM   1166  O3*   U    16      29.263  -0.099   5.331  1.00  0.00           O  
ATOM   1167  C2*   U    16      27.419   1.134   4.239  1.00  0.00           C  
ATOM   1168  O2*   U    16      26.823   1.075   5.541  1.00  0.00           O  
ATOM   1169  C1*   U    16      27.333   2.547   3.658  1.00  0.00           C  
ATOM   1170  N1    U    16      27.040   2.530   2.195  1.00  0.00           N  
ATOM   1171  C2    U    16      25.709   2.427   1.812  1.00  0.00           C  
ATOM   1172  O2    U    16      24.794   2.353   2.631  1.00  0.00           O  
ATOM   1173  N3    U    16      25.464   2.408   0.449  1.00  0.00           N  
ATOM   1174  C4    U    16      26.416   2.482  -0.553  1.00  0.00           C  
ATOM   1175  O4    U    16      26.084   2.456  -1.738  1.00  0.00           O  
ATOM   1176  C5    U    16      27.775   2.588  -0.069  1.00  0.00           C  
ATOM   1177  C6    U    16      28.043   2.610   1.260  1.00  0.00           C  
ATOM   1178 1H5*   U    16      31.393   3.204   4.777  1.00  0.00           H  
ATOM   1179 2H5*   U    16      31.482   1.464   4.432  1.00  0.00           H  
ATOM   1180  H4*   U    16      29.391   2.412   5.717  1.00  0.00           H  
ATOM   1181  H3*   U    16      29.282   0.519   3.323  1.00  0.00           H  
ATOM   1182 1H2*   U    16      26.932   0.424   3.564  1.00  0.00           H  
ATOM   1183 2HO*   U    16      26.287   0.279   5.574  1.00  0.00           H  
ATOM   1184  H1*   U    16      26.554   3.105   4.182  1.00  0.00           H  
ATOM   1185  H3    U    16      24.499   2.333   0.158  1.00  0.00           H  
ATOM   1186  H5    U    16      28.596   2.645  -0.783  1.00  0.00           H  
ATOM   1187  H6    U    16      29.078   2.699   1.589  1.00  0.00           H  
ATOM   1188  P     A    17      29.128  -1.691   5.090  1.00  0.00           P  
ATOM   1189  O1P   A    17      30.089  -2.375   5.985  1.00  0.00           O  
ATOM   1190  O2P   A    17      29.170  -1.952   3.633  1.00  0.00           O  
ATOM   1191  O5*   A    17      27.632  -2.003   5.614  1.00  0.00           O  
ATOM   1192  C5*   A    17      27.278  -1.850   6.998  1.00  0.00           C  
ATOM   1193  C4*   A    17      25.805  -2.202   7.251  1.00  0.00           C  
ATOM   1194  O4*   A    17      24.933  -1.552   6.313  1.00  0.00           O  
ATOM   1195  C3*   A    17      25.537  -3.691   7.076  1.00  0.00           C  
ATOM   1196  O3*   A    17      25.848  -4.377   8.302  1.00  0.00           O  
ATOM   1197  C2*   A    17      24.052  -3.687   6.737  1.00  0.00           C  
ATOM   1198  O2*   A    17      23.249  -3.687   7.926  1.00  0.00           O  
ATOM   1199  C1*   A    17      23.841  -2.405   5.932  1.00  0.00           C  
ATOM   1200  N9    A    17      23.844  -2.667   4.467  1.00  0.00           N  
ATOM   1201  C8    A    17      24.786  -2.368   3.532  1.00  0.00           C  
ATOM   1202  N7    A    17      24.548  -2.718   2.311  1.00  0.00           N  
ATOM   1203  C5    A    17      23.294  -3.330   2.427  1.00  0.00           C  
ATOM   1204  C6    A    17      22.444  -3.936   1.494  1.00  0.00           C  
ATOM   1205  N6    A    17      22.737  -4.034   0.196  1.00  0.00           N  
ATOM   1206  N1    A    17      21.280  -4.436   1.948  1.00  0.00           N  
ATOM   1207  C2    A    17      20.966  -4.348   3.242  1.00  0.00           C  
ATOM   1208  N3    A    17      21.695  -3.799   4.210  1.00  0.00           N  
ATOM   1209  C4    A    17      22.857  -3.304   3.734  1.00  0.00           C  
ATOM   1210 1H5*   A    17      27.451  -0.815   7.297  1.00  0.00           H  
ATOM   1211 2H5*   A    17      27.911  -2.502   7.603  1.00  0.00           H  
ATOM   1212  H4*   A    17      25.528  -1.907   8.258  1.00  0.00           H  
ATOM   1213  H3*   A    17      26.119  -4.088   6.237  1.00  0.00           H  
ATOM   1214 1H2*   A    17      23.812  -4.551   6.120  1.00  0.00           H  
ATOM   1215 2HO*   A    17      22.338  -3.821   7.656  1.00  0.00           H  
ATOM   1216  H1*   A    17      22.896  -1.943   6.222  1.00  0.00           H  
ATOM   1217  H8    A    17      25.704  -1.843   3.798  1.00  0.00           H  
ATOM   1218 1H6    A    17      22.086  -4.481  -0.435  1.00  0.00           H  
ATOM   1219 2H6    A    17      23.607  -3.662  -0.156  1.00  0.00           H  
ATOM   1220  H2    A    17      20.003  -4.767   3.538  1.00  0.00           H  
ATOM   1221  P     U    18      26.579  -5.817   8.301  1.00  0.00           P  
ATOM   1222  O1P   U    18      26.793  -6.227   9.708  1.00  0.00           O  
ATOM   1223  O2P   U    18      27.729  -5.760   7.371  1.00  0.00           O  
ATOM   1224  O5*   U    18      25.456  -6.782   7.657  1.00  0.00           O  
ATOM   1225  C5*   U    18      24.254  -7.102   8.372  1.00  0.00           C  
ATOM   1226  C4*   U    18      23.146  -7.595   7.431  1.00  0.00           C  
ATOM   1227  O4*   U    18      22.934  -6.690   6.339  1.00  0.00           O  
ATOM   1228  C3*   U    18      23.513  -8.921   6.759  1.00  0.00           C  
ATOM   1229  O3*   U    18      23.167 -10.054   7.578  1.00  0.00           O  
ATOM   1230  C2*   U    18      22.670  -8.867   5.493  1.00  0.00           C  
ATOM   1231  O2*   U    18      21.344  -9.360   5.732  1.00  0.00           O  
ATOM   1232  C1*   U    18      22.649  -7.384   5.114  1.00  0.00           C  
ATOM   1233  N1    U    18      23.642  -7.059   4.047  1.00  0.00           N  
ATOM   1234  C2    U    18      23.311  -7.401   2.742  1.00  0.00           C  
ATOM   1235  O2    U    18      22.244  -7.945   2.458  1.00  0.00           O  
ATOM   1236  N3    U    18      24.246  -7.089   1.771  1.00  0.00           N  
ATOM   1237  C4    U    18      25.468  -6.473   1.981  1.00  0.00           C  
ATOM   1238  O4    U    18      26.225  -6.242   1.039  1.00  0.00           O  
ATOM   1239  C5    U    18      25.737  -6.149   3.366  1.00  0.00           C  
ATOM   1240  C6    U    18      24.838  -6.444   4.339  1.00  0.00           C  
ATOM   1241 1H5*   U    18      23.903  -6.212   8.897  1.00  0.00           H  
ATOM   1242 2H5*   U    18      24.473  -7.881   9.105  1.00  0.00           H  
ATOM   1243  H4*   U    18      22.215  -7.704   7.990  1.00  0.00           H  
ATOM   1244  H3*   U    18      24.579  -8.932   6.500  1.00  0.00           H  
ATOM   1245 1H2*   U    18      23.153  -9.447   4.704  1.00  0.00           H  
ATOM   1246 2HO*   U    18      21.407 -10.312   5.833  1.00  0.00           H  
ATOM   1247  H1*   U    18      21.646  -7.113   4.775  1.00  0.00           H  
ATOM   1248  H3    U    18      24.014  -7.334   0.818  1.00  0.00           H  
ATOM   1249  H5    U    18      26.676  -5.662   3.631  1.00  0.00           H  
ATOM   1250  H6    U    18      25.068  -6.171   5.372  1.00  0.00           H  
ATOM   1251  P     G    19      24.070 -11.395   7.584  1.00  0.00           P  
ATOM   1252  O1P   G    19      23.816 -12.113   8.853  1.00  0.00           O  
ATOM   1253  O2P   G    19      25.456 -11.030   7.213  1.00  0.00           O  
ATOM   1254  O5*   G    19      23.435 -12.255   6.371  1.00  0.00           O  
ATOM   1255  C5*   G    19      22.522 -13.339   6.622  1.00  0.00           C  
ATOM   1256  C4*   G    19      21.910 -13.906   5.338  1.00  0.00           C  
ATOM   1257  O4*   G    19      21.631 -12.876   4.378  1.00  0.00           O  
ATOM   1258  C3*   G    19      22.872 -14.848   4.607  1.00  0.00           C  
ATOM   1259  O3*   G    19      22.849 -16.193   5.123  1.00  0.00           O  
ATOM   1260  C2*   G    19      22.306 -14.797   3.194  1.00  0.00           C  
ATOM   1261  O2*   G    19      21.217 -15.718   3.035  1.00  0.00           O  
ATOM   1262  C1*   G    19      21.836 -13.351   3.036  1.00  0.00           C  
ATOM   1263  N9    G    19      22.839 -12.528   2.309  1.00  0.00           N  
ATOM   1264  C8    G    19      23.791 -11.684   2.790  1.00  0.00           C  
ATOM   1265  N7    G    19      24.552 -11.088   1.933  1.00  0.00           N  
ATOM   1266  C5    G    19      24.065 -11.583   0.719  1.00  0.00           C  
ATOM   1267  C6    G    19      24.486 -11.309  -0.611  1.00  0.00           C  
ATOM   1268  O6    G    19      25.390 -10.564  -0.986  1.00  0.00           O  
ATOM   1269  N1    G    19      23.731 -12.014  -1.541  1.00  0.00           N  
ATOM   1270  C2    G    19      22.696 -12.877  -1.236  1.00  0.00           C  
ATOM   1271  N2    G    19      22.087 -13.463  -2.269  1.00  0.00           N  
ATOM   1272  N3    G    19      22.295 -13.140   0.012  1.00  0.00           N  
ATOM   1273  C4    G    19      23.016 -12.464   0.937  1.00  0.00           C  
ATOM   1274 1H5*   G    19      21.699 -12.982   7.229  1.00  0.00           H  
ATOM   1275 2H5*   G    19      23.061 -14.131   7.167  1.00  0.00           H  
ATOM   1276  H4*   G    19      20.985 -14.435   5.577  1.00  0.00           H  
ATOM   1277  H3*   G    19      23.888 -14.436   4.620  1.00  0.00           H  
ATOM   1278 1H2*   G    19      23.095 -15.012   2.471  1.00  0.00           H  
ATOM   1279 2HO*   G    19      20.988 -15.731   2.103  1.00  0.00           H  
ATOM   1280  H1*   G    19      20.886 -13.332   2.498  1.00  0.00           H  
ATOM   1281  H8    G    19      23.908 -11.512   3.860  1.00  0.00           H  
ATOM   1282  H1    G    19      23.976 -11.870  -2.509  1.00  0.00           H  
ATOM   1283 1H2    G    19      21.324 -14.104  -2.107  1.00  0.00           H  
ATOM   1284 2H2    G    19      22.390 -13.265  -3.212  1.00  0.00           H  
ATOM   1285  P     C    20      24.094 -17.201   4.902  1.00  0.00           P  
ATOM   1286  O1P   C    20      23.887 -18.383   5.769  1.00  0.00           O  
ATOM   1287  O2P   C    20      25.351 -16.425   5.009  1.00  0.00           O  
ATOM   1288  O5*   C    20      23.916 -17.661   3.361  1.00  0.00           O  
ATOM   1289  C5*   C    20      23.025 -18.728   3.003  1.00  0.00           C  
ATOM   1290  C4*   C    20      22.893 -18.887   1.479  1.00  0.00           C  
ATOM   1291  O4*   C    20      22.785 -17.619   0.811  1.00  0.00           O  
ATOM   1292  C3*   C    20      24.132 -19.534   0.852  1.00  0.00           C  
ATOM   1293  O3*   C    20      24.116 -20.969   0.981  1.00  0.00           O  
ATOM   1294  C2*   C    20      23.970 -19.088  -0.597  1.00  0.00           C  
ATOM   1295  O2*   C    20      23.083 -19.960  -1.314  1.00  0.00           O  
ATOM   1296  C1*   C    20      23.395 -17.672  -0.492  1.00  0.00           C  
ATOM   1297  N1    C    20      24.448 -16.628  -0.664  1.00  0.00           N  
ATOM   1298  C2    C    20      24.788 -16.267  -1.967  1.00  0.00           C  
ATOM   1299  O2    C    20      24.233 -16.804  -2.926  1.00  0.00           O  
ATOM   1300  N3    C    20      25.743 -15.314  -2.147  1.00  0.00           N  
ATOM   1301  C4    C    20      26.346 -14.734  -1.103  1.00  0.00           C  
ATOM   1302  N4    C    20      27.278 -13.806  -1.323  1.00  0.00           N  
ATOM   1303  C5    C    20      26.005 -15.096   0.239  1.00  0.00           C  
ATOM   1304  C6    C    20      25.058 -16.040   0.414  1.00  0.00           C  
ATOM   1305 1H5*   C    20      22.039 -18.520   3.422  1.00  0.00           H  
ATOM   1306 2H5*   C    20      23.399 -19.662   3.428  1.00  0.00           H  
ATOM   1307  H4*   C    20      22.009 -19.486   1.252  1.00  0.00           H  
ATOM   1308  H3*   C    20      25.046 -19.104   1.281  1.00  0.00           H  
ATOM   1309 1H2*   C    20      24.948 -19.058  -1.085  1.00  0.00           H  
ATOM   1310 2HO*   C    20      22.206 -19.849  -0.939  1.00  0.00           H  
ATOM   1311  H1*   C    20      22.621 -17.539  -1.252  1.00  0.00           H  
ATOM   1312 1H4    C    20      27.743 -13.361  -0.545  1.00  0.00           H  
ATOM   1313 2H4    C    20      27.521 -13.549  -2.269  1.00  0.00           H  
ATOM   1314  H5    C    20      26.495 -14.624   1.092  1.00  0.00           H  
ATOM   1315  H6    C    20      24.775 -16.334   1.424  1.00  0.00           H  
ATOM   1316  P     C    21      25.478 -21.837   0.913  1.00  0.00           P  
ATOM   1317  O1P   C    21      25.166 -23.216   1.355  1.00  0.00           O  
ATOM   1318  O2P   C    21      26.555 -21.078   1.590  1.00  0.00           O  
ATOM   1319  O5*   C    21      25.806 -21.872  -0.669  1.00  0.00           O  
ATOM   1320  C5*   C    21      25.064 -22.713  -1.566  1.00  0.00           C  
ATOM   1321  C4*   C    21      25.331 -22.359  -3.037  1.00  0.00           C  
ATOM   1322  O4*   C    21      25.492 -20.944  -3.228  1.00  0.00           O  
ATOM   1323  C3*   C    21      26.646 -22.955  -3.545  1.00  0.00           C  
ATOM   1324  O3*   C    21      26.476 -24.311  -3.980  1.00  0.00           O  
ATOM   1325  C2*   C    21      26.966 -22.031  -4.712  1.00  0.00           C  
ATOM   1326  O2*   C    21      26.289 -22.448  -5.905  1.00  0.00           O  
ATOM   1327  C1*   C    21      26.470 -20.658  -4.244  1.00  0.00           C  
ATOM   1328  N1    C    21      27.584 -19.805  -3.730  1.00  0.00           N  
ATOM   1329  C2    C    21      28.285 -19.039  -4.659  1.00  0.00           C  
ATOM   1330  O2    C    21      27.988 -19.086  -5.851  1.00  0.00           O  
ATOM   1331  N3    C    21      29.300 -18.247  -4.213  1.00  0.00           N  
ATOM   1332  C4    C    21      29.619 -18.201  -2.914  1.00  0.00           C  
ATOM   1333  N4    C    21      30.618 -17.412  -2.517  1.00  0.00           N  
ATOM   1334  C5    C    21      28.906 -18.983  -1.951  1.00  0.00           C  
ATOM   1335  C6    C    21      27.903 -19.767  -2.396  1.00  0.00           C  
ATOM   1336 1H5*   C    21      23.999 -22.597  -1.362  1.00  0.00           H  
ATOM   1337 2H5*   C    21      25.345 -23.754  -1.392  1.00  0.00           H  
ATOM   1338  H4*   C    21      24.502 -22.712  -3.654  1.00  0.00           H  
ATOM   1339  H3*   C    21      27.425 -22.882  -2.777  1.00  0.00           H  
ATOM   1340 1H2*   C    21      28.048 -22.003  -4.876  1.00  0.00           H  
ATOM   1341 2HO*   C    21      26.713 -23.255  -6.208  1.00  0.00           H  
ATOM   1342  H1*   C    21      25.977 -20.150  -5.077  1.00  0.00           H  
ATOM   1343 1H4    C    21      30.870 -17.368  -1.539  1.00  0.00           H  
ATOM   1344 2H4    C    21      31.124 -16.859  -3.193  1.00  0.00           H  
ATOM   1345  H5    C    21      29.165 -18.946  -0.892  1.00  0.00           H  
ATOM   1346  H6    C    21      27.342 -20.373  -1.686  1.00  0.00           H  
ATOM   1347  H3T   C    21      27.352 -24.670  -4.145  1.00  0.00           H  
TER    1348        C    21                                                      
ENDMDL                                                                  
"""


if __name__ == '__main__':
    main()

