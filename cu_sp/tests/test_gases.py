import sys
import os
sys.path.insert(0, os.path.abspath(r'..'))

import gases
import unittest

props, units = gases.get_properties()
water = props.loc['H2O']

class TestSequenceFunctions(unittest.TestCase):
    def test_H2O(self):
        self.assertTrue(water.TriplePointT == 273.15)
    
    def test_units(self):
        self.assertTrue( units['L_fusion'] == 'J/kg')
        
if __name__ == '__main__':
    unittest.main()
