# -*- coding: utf-8 -*-
import rlcompleter, readline
readline.parse_and_bind('tab:complete')
#############################""
import unittest
from Tools import *
import numpy as np


class GuillaumeTest(unittest.TestCase):  
  def setUp(self):
    pass  

  def tearDown(self):
    pass

  def testXYZ(self):
    a = 1
    b = 2
    f = True
    self.assertEqual(3, a+b)
    self.assertTrue(f)
    #self.assertRaises(ValueError, self.toto, 1, 2)
    return
    
  def testABC(self):
      return
      
  def assertNumpyEqual(self, a, b):
      diff = np.sum(np.abs(a-b)) 
      return diff < 1.0e-6
  
  
  def testScalarFieldFromFile(self):
      x = np.arange(0,1,0.1)
      y = x*x
      np.savetxt("test.dt_ev", np.c_[x,y], header="x y")
      f = Field.LoadFromFile("test.dt_ev", ["y"])
      self.assertTrue(self.assertNumpyEqual(y, f._npa))
      return

if __name__ == "__main__":
  suite = unittest.TestSuite()
#   suite.addTest(GuillaumeTest('testXYZ'))
#   suite.addTest(GuillaumeTest('testABC'))
  suite.addTest(GuillaumeTest('testScalarFieldFromFile'))
  unittest.TextTestRunner().run(suite)


