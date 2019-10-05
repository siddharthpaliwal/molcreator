"""
Calls the writer modules to
generate the coordinates and write a LAMMPS input file
"""

from molcreator.system import System
from molcreator.molecule import Molecule

if __name__ == '__main__':
    path = './test_output/'
    boxlen = [20.0, 20.0, 20.0]
    normal = [0.0, 0.0, 1.0]
    decane = Molecule('DEC', natoms=10, nbonds=9, nangles=8, ndihedrals=7)
    system = System(decane, nmol=100, box=boxlen)
    system.gen_manifold('planar', normal=[0.0, 0.0, 1.0], tol=1.6)
    system.gen_molecules(decane)
    print('System Box', system.box)
    system.write_settings(path)
    system.write_coords_lmp(path)
