# -*- coding: utf-8 -*-
"""
Created on Tue Nov  4 19:01:54 2014

@author: Dmitri
"""
import unittest
from population import *
import numpy as np


class TestPopulation(unittest.TestCase):
    def testStationary(self):
        MSG = "stationary vector must be in the null-space of the matrix Q!"
        DELTA = 1e-12
        for N in range(2, 30):
            pop = population(N);
            p_Q = np.linalg.norm(np.dot(pop.Pstat, pop.Q))
            self.assertAlmostEqual(p_Q, 0, places=None, msg=MSG, delta=DELTA)

    def testNegative(self):
        """should fail with negative input"""
        self.assertRaises(OutOfRangeError, population, -1)


if __name__ == "__main__":
    unittest.main() 