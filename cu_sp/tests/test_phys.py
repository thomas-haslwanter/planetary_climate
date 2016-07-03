import sys
import os
sys.path.insert(0, os.path.abspath(r'..'))

import phys
import unittest
import numpy as np
from scipy import constants

class TestSequenceFunctions(unittest.TestCase):
    def test_const(self):
        self.assertAlmostEqual( phys.c - constants.c, 0 )
    
    def test_Plank(self):
        self.assertAlmostEqual( phys.B(1.5e13, 300), 4.966991e-12 )
        
if __name__ == '__main__':
    unittest.main()
