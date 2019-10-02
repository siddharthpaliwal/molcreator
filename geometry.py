from molecule import Molecule
import numpy as np
import scipy.spatial.distance as spdist
import time


class Geometry:
    """
    Define of some common geometrical manifolds for arranging molecules
    """
    box = np.zeros((3, 2))  # [ [xlo,xhi], [ylo,yhi], [zlo,zhi] ]

    def __init__(self,
                 origin: (float, float, float) = (0.0, 0.0, 0.0),
                 boxlen: (float, float, float) = (1.0, 1.0, 1.0)):
        """
        Initialize a geometric manifold at specified origin
        :param origin: Origin of the geometric manifold, default (0, 0, 0)
        :param box: Box length in each dimension, default (1.0, 1.0, 1.0)
        """
        self.origin = origin
        self.boxlen = boxlen


class Planar(Geometry):
    """
    Define a planar manifold
    """

    def __init__(self,
                 npts,
                 origin: (float, float, float) = (0.0, 0.0, 0.0),
                 boxlen: (float, float, float) = (1.0, 1.0, 1.0),
                 normal: (float, float, float) = (0.0, 0.0, 1.0)):
        """
        Define a planar manifold containing npts points
        :param npts: No. of seed points
        :param origin: Box origin
        :param boxlen: Boxlength in each dim. default:(1.0, 1.0, 1.0)
        :param normal: Vector normal to the planar manifold
        """
        super().__init__(origin, boxlen)
        self.boxlen[2] = 0.0  # make z dim = 0
        self.normal = normal
        self.seeds = np.array([[0.0, 0.0, 0.0] for _ in range(npts)])
        # self.seeds = [[0.0, 0.0] for _ in range(npts)]

    def generate_seeds(self, tol=1.0):
        """
        Generate seed points on the manifold atleast 'tol' distance apart
        """
        _nassigned = 1  # assigned seeds
        self.seeds[0] = np.random.rand(3) * self.boxlen
        _npts = len(self.seeds)
        _ntry = 0  # trial attempts for each seed
        start = time.time()
        while _nassigned < _npts:
            rtest = np.random.rand(3) * self.boxlen
            dist = spdist.cdist(self.seeds[:_nassigned], np.asarray([rtest]))
            if _ntry > 1e2:
                print("Failed at too many attempts, try changing tol value or box size")
                break
            if np.any(dist) < tol:
                _ntry += 1
                continue
            else:
                self.seeds[_nassigned] = rtest
                _nassigned += 1
                _ntry = 0
                print("Assigned {:d}".format(_nassigned))
            # trial_seeds = np.random.uniform((nrem,2))*np.tile(self.boxlen,(nrem,1))
            # distmat = spdist.squareform(spdist.pdist(trial_seeds))
        end = time.time()
        print('Generated {:d} seeds in {:.3f} sec.'.format(_nassigned, (end - start)))
        # [print(seed) for seed in self.seeds]


if __name__ == '__main__':
    plane = Planar(10000, boxlen=[5.0, 5.0, 5.0])
    plane.generate_seeds(0.1)
