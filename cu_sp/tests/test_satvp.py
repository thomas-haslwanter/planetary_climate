import sys
import os
sys.path.insert(0, os.path.abspath(r'..'))

import satvp
import unittest
import numpy as np

class TestSequenceFunctions(unittest.TestCase):
    def test_H2O(self):
        self.assertAlmostEqual(satvp.satvp_H2O(300), 3589.9143379302436)
        self.assertAlmostEqual(satvp.satvp_H2O(260, mode='ice'), 195.4964678727905)
        self.assertEqual(satvp.satvp_H2O(260),  satvp.satvp_H2O(260, mode='general') )
        self.assertEqual( len(satvp.satvp_H2O(np.arange(250,260, 2))), 5)
    '''
    [To be done]
    def test_satvp(self):
    
    def test_MoistAdiabat(self):
    '''
        
if __name__ == '__main__':
    unittest.main()
