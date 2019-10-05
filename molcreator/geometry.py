import numpy as np
import scipy.spatial.distance as spdist
import time
import copy
from molcreator.poisson_sampling import PoissonDisc


class Geometry:
    """Define of some common geometrical manifolds for arranging molecules
    """

    def __init__(self,
                 origin: [float, float, float] = (0.0, 0.0, 0.0),
                 boxlen: [float, float, float] = (1.0, 1.0, 1.0)):
        """Initialize a geometric manifold at specified origin

        Args:
            origin: Origin of the geometric manifold, default (0, 0, 0)
            boxlen: Box length in each dimension, default (1.0, 1.0, 1.0)
        """
        self.origin = copy.copy(origin)
        self.boxlen = copy.copy(boxlen)


class Planar(Geometry):
    """Define a planar manifold

    Args:
        npts: No. of seed points
        origin: Box origin
        boxlen: Boxlength in each dim. default:(1.0, 1.0, 1.0)
        normal (float, float, float): Vector normal to the planar manifold

    Attributes:
        origin: Box origin
        boxlen: Boxlength in each dim. default:(1.0, 1.0, 1.0)
        normal: Vector normal to the planar manifold
    """

    def __init__(self,
                 npts,
                 origin: [float, float, float] = (0.0, 0.0, 0.0),
                 boxlen: [float, float, float] = (1.0, 1.0, 1.0),
                 normal: [float, float, float] = (0.0, 0.0, 1.0)):
        """Define a planar manifold containing npts points

        Args:
            npts: No. of seed points
            origin: Box origin
            boxlen: Boxlength in each dim. default:(1.0, 1.0, 1.0)
            normal (float, float, float): Vector normal to the planar manifold
        """
        super().__init__(origin, boxlen)
        self.boxlen[2] = 0.0  # make z dim = 0
        self.normal = normal
        self.seeds = np.array([[0.0, 0.0, 0.0] for _ in range(npts)])
        # self.seeds = [[0.0, 0.0] for _ in range(npts)]

    def generate_seeds(self, tol=1.0):
        """Generate seed points on the manifold atleast 'tol' distance apart

        Args:
            tol: Min. distance between any two seeds
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

    def generate_seeds_poisson(self, tol=1.0):
        npts = len(self.seeds)
        pdisc = PoissonDisc(self.boxlen[0], self.boxlen[1], r=tol)
        start = time.time()
        coords = np.asarray(pdisc.sample())
        end = time.time()
        print(np.shape(coords))
        print('Generated {:d} seeds in {:.3f} sec.'.format(len(coords), (end - start)))
        if len(coords) >= npts:
            pts = np.random.randint(0, len(coords), npts)
            self.seeds = np.column_stack((coords[pts], np.zeros(npts)))
        else:
            print('Couldnt generate {:d} coords. Reduce tol value'.format(npts))
            self.seeds = np.zeros((npts, 3))

    def generate_seeds_square(self, tol=1.0):
        npts = len(self.seeds)
        start = time.time()
        nx = ny = np.ceil(np.sqrt(npts))
        xlist = ylist = np.linspace(self.origin[0], self.boxlen[0], nx, endpoint=False)
        xcoords, ycoords = np.meshgrid(xlist, ylist)
        end = time.time()
        ngen = np.size(xcoords)
        # print(np.shape(xcoords), xcoords[:10])
        print('Generated {:d} seeds in {:.3f} sec.'.format(np.size(xcoords), (end - start)))
        if ngen > npts:
            pts = np.random.randint(0, ngen, npts)
            self.seeds = np.column_stack((np.reshape(xcoords,(ngen, 1))[pts], np.reshape(ycoords,(ngen, 1))[pts], np.zeros(npts)))
        elif ngen == npts:
            self.seeds = np.column_stack((np.reshape(xcoords, (ngen, 1)), np.reshape(ycoords, (ngen, 1)), np.zeros(npts)))
        else:
            print('Couldnt generate {:d} coords. Reduce tol value'.format(npts))
            self.seeds = np.zeros((npts, 3))


if __name__ == '__main__':
    plane = Planar(100, boxlen=[5.0, 5.0, 5.0])
    # plane.generate_seeds(0.1)
    plane.generate_seeds_poisson(0.4)
