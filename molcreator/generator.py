"""
Calls the writer modules to
generate the coordinates and write a LAMMPS input file
"""

from molcreator.system import System, write_settings

if __name__ == '__main__':
    path = './test_output/'
    boxlen = [20.0, 20.0, 20.0]
    normal = [0.0, 0.0, 1.0]
    system = System(nmol=10, box=boxlen, geom_type='planar',
                    moltype='DEC', natompermol=10)
    print("Now generating seeds and molecules")
    write_settings(path)
    system.write_coords_lmp(path)
