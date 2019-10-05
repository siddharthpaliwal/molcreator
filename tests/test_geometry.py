import unittest
from molcreator.geometry import Planar
import numpy as np


class TestGeometry(unittest.TestCase):
    def test_init(self):
        testbox = [5.0, 5.0, 5.0]
        plane = Planar(100, boxlen=testbox)
        self.assertListEqual(plane.boxlen[:2], testbox[:2], "Box not assigned correctly")

    def test_normal(self):
        testnormal = [0.5, 0.5, 0.0]
        plane = Planar(100, normal=testnormal)
        self.assertListEqual(plane.normal, testnormal)

    def test_generate_seed(self):
        tol = 0.1
        plane = Planar(100, boxlen=[5.0, 5.0, 5.0])
        plane.generate_seeds(tol)
        randi = np.random.randint(0, 100, 2)
        dist = np.sum((plane.seeds[randi[0]] - plane.seeds[randi[1]]) ** 2) ** 0.5
        self.assertLessEqual(tol, dist, "distance < tolerance")


if __name__ == '__main__':
    unittest.main()
