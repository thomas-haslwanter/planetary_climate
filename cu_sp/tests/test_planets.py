import sys
import os
sys.path.insert(0, os.path.abspath(r'..'))

import planets
import unittest

props, units = planets.get_properties()
earth = props.loc['Earth']

class TestSequenceFunctions(unittest.TestCase):
    def test_earth(self):
        self.assertTrue(earth.year == 365.256 * 3600 * 24)
    
    def test_units(self):
        self.assertTrue( units['year'] == 'sec')
        
if __name__ == '__main__':
    unittest.main()
