'''
Unit tests for exported function of ace_cream.py
'''

import unittest
import numpy as np
from ace_cream import f_mapping
class Test_Ace_Cream(unittest.TestCase):
    def test_f_mapping(self):
        x = [1,2,1,3]
        tx = np.array([[1,1],[2,2],[1,1],[3,3]])
        x_candidate=[2,1]
        t_x_candidate = np.array([[2,2],[1,1]])
        return (t_x_candidate == f_mapping(x, tx, x_candidate)).all()
      
if __name__ == "__main__":
    unittest.main()