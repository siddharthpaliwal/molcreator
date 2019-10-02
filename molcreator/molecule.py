"""
The :class:`~.Molecule` defines the basic unit of a molecule in the system.
The attributes are defined to correspond to a United Atom force field model
for a linear chain alkane.
"""

import numpy as np
from molcreator.trappe import Trappe
from math import sin, cos, pi
from scipy.spatial.transform import Rotation


class Molecule:
    R"""Defines the molecule class used to instantiate a molecule objects
    that could be modified in space and arranged in the desired manner

    Args:
        name (string): Chemical/Ref. name of the molecule
        natoms (integer): No. of beads in the molecule representation

    Attributes:
        length (float): Max. length of the molecule
        coords ()
    """
    name = ''
    natoms: int = 0
    beadtypes = ('CH4', 'CH3', 'CH2')
    type = []
    bond = []

    def __init__(self, name: str, natoms: int):
        """

        :type name: string
        :type natoms: integer
        """
        self.name = name
        self.natoms = natoms
        self.paxis = (1.0, 0.0, 0.0)  # principal axis
        self.base_coords = (0.0, 0.0, 0.0)  # base point of the elongated molecule
        self.atom_coords = np.zeros((natoms, 3))
        self.nbonds = None

        # Call the initialization functions
        self.gen_coords()
        self.gen_types()
        self.gen_bonds()

    def gen_coords(self):
        """
        Generate the coordinates of atoms in the molecule
        :return:
        """
        bondlength = Trappe['bondlength']
        bondangle = Trappe['angles']['backbone']['theta0'] * pi / 180.
        # Assign coordinates
        for atom in range(0, self.natoms, 2):
            self.atom_coords[atom, 0] = atom * bondlength * sin(bondangle / 2)
            self.atom_coords[atom, 1] = 0.0
            self.atom_coords[atom, 2] = 0.0
        for atom in range(1, self.natoms, 2):
            self.atom_coords[atom, 0] = atom * bondlength * sin(bondangle / 2)
            self.atom_coords[atom, 1] = cos(bondangle / 2)
            self.atom_coords[atom, 2] = 0.0

        return

    def gen_types(self):
        """
        Assign the types to each atom in the molecule
        """
        self.type.append('CH3')
        for atom in range(1, self.natoms):
            self.type.append('CH2')
        self.type.append('CH3')

    def gen_bonds(self):
        """
        Assign the bonds and bond-type to the molecule
        """
        for atomi, atomj in zip(range(1,self.natoms-1), range(2,self.natoms)):
            self.bond.append(('sat', atomi, atomj))
        return

    def print_coords(self):
        """
        Print the coordinates of atoms in the molecule for debugging
        :return: None
        """
        for atom in range(self.natoms):
            print('Atom:{:d} Type:{:s} {:.3f} {:.3f} {:.3f}'.format(atom, self.type[atom],
                                                                    self.atom_coords[atom, 0],
                                                                    self.atom_coords[atom, 1],
                                                                    self.atom_coords[atom, 2]))

    def translate(self, delta: (float, float, float)):
        """
        Translate the molecule by given distance in x,y,z directions
        :return: A copy of the translated coordinates
        """
        self.atom_coords += np.tile(delta, (self.natoms, 1))
        return self.atom_coords + np.tile(delta, (self.natoms, 1))
    
    def move_to_coords(self, pos: (float, float, float)) -> None:
        """
        Move the molecule to the given position
        """
        self.translate(pos-self.base_coords)
        self.base_coords = pos
        return

    # def check_overlap(self, rtest: (float, float, float)):
    #     for i in range(self.nmol):
    #         dr = np.sum((rtest - self.atom_coords[i])**2)**0.5
    #         if(dr<self.tolerance):
    #             return True
    #     return False

    def get_distance(self, moltest):
        """
        Calculate the min. distance between any two atoms of self from moltest
        """
        return

    def check_overlap(self, ri, rj):
        """
        Check if the two molecules overlap
        """
        return True

    def rotate_paxis(self, axis) -> None:
        """
        rotate the molecule such that the principle axis is now pointing in 'axis' direction
        """
        nvec = np.cross(self.paxis, axis) # normal vector to current and target axis
        theta = np.arcsin(np.linalg.norm(nvec)) # angle between current and target axis
        norm = np.linalg.norm(nvec)
        if(norm > 0.0):
            nvec = nvec/norm # normalize
        r = Rotation.from_rotvec(theta*nvec) # get the rotation matrix
        self.atom_coords = r.apply(self.atom_coords) # apply the rotation to each atom of the molecule
        self.paxis = axis
        return

