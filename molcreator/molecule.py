"""
The :class:`~.Molecule` defines the basic unit of a molecule in the system.
The attributes are defined to correspond to a United Atom force field model
for a linear chain alkane.
"""

import numpy as np
from math import sin, cos, pi
from scipy.spatial.transform import Rotation
from molcreator.trappe import Trappe, AtomType, BondType, AngleType, DihedralType


class Atom:
    """Stores the information of each atom

    Attributes:
        type (:class AtomType): Type of Atom
        coords (float, float, float): coordinates of the atom
        index (int): index
    """

    def __init__(self, coords=[0.0, 0.0, 0.0], atype=AtomType.CH2):
        self.type = AtomType(atype)
        self.coords = coords
        self.index = 0


class Bond:
    """Stores the information of each bond

    Attributes:
        type (:class BondType): Type of Bond
        atoms (int, int): tuple of 2 consecutive atoms connected by the bond
        index (int): index
    """

    def __init__(self, atomi: int, atomj: int, btype=BondType.SAT):
        """

        Args:
            btype (object): BondType Enum
        """
        self.type = BondType(btype)
        self.atoms = (atomi, atomj)
        self.index = 0


class Angle:
    """Stores the information of each angle

    Attributes:
        type: Type of Angle
        atoms (int, int, int): tuple of 3 consecutive atoms forming the angle
        index (int): index
    """

    def __init__(self, atomi: int, atomj: int, atomk: int, atype=AngleType.BACKBONE):
        self.type = AngleType(atype)
        self.atoms = (atomi, atomj, atomk)
        self.index = 0


class Dihedral:
    """Stores the information of each dihedral

    Attributes:
        type: Type of Dihedral angle
        atoms (int, int, int, int): tuple of 4 consecutive atoms forming the dihedral angle
        index (int): index
    """

    def __init__(self, atomi: int, atomj: int, atomk: int, atoml: int, atype=DihedralType.BACKBONE):
        self.type = DihedralType(atype)
        self.atoms = (atomi, atomj, atomk, atoml)
        self.index = 0


class Molecule:
    R"""Defines the molecule class used to instantiate a molecule objects
    that could be modified in space and arranged in the desired manner

    Args:
        name (string): Chemical/Ref. name of the molecule
        natoms (integer): No. of beads in the molecule representation
        nbonds (integer): No. of beads in the molecule representation
        nangles (integer): No. of beads in the molecule representation
        ndihedrals (integer): No. of beads in the molecule representation

    Attributes:
        base_coords: coordinates of the base
        paxis: principle axis unit vector in 3D
        index: index number of the molecule
        atoms (list): list of :class: Atom objects
        bonds (list): list of :class: Bond objects
        angles (list): list of :class: Angle objects
        dihedrals (list): list of :class: Dihedral objects
    """

    def __init__(self,
                 name: str,
                 natoms: int,
                 nbonds: int,
                 nangles: int,
                 ndihedrals: int,
                 idx: int = 0):
        """
        Initialize a molecule instance
        Args:
            name: 3 letter name of the molecule
            natoms: NO. of atoms in the molecule
        """
        assert (len(name) == 3 and isinstance(name, str))
        self.name = name
        self.natoms = natoms
        self.nbonds = nbonds
        self.nangles = nangles
        self.ndihedrals = ndihedrals
        self.paxis = (1.0, 0.0, 0.0)  # principal axis
        self.base_coords = (0.0, 0.0, 0.0)  # base point of the elongated molecule
        self.index = idx
        # self.atom_coords = np.zeros((natoms, 3))
        # self.atoms = [{'type':'', 'coords':[0.0, 0.0, 0.0]} for _ in range(natoms)]
        self.atoms = [Atom() for _ in range(self.natoms)]
        self.bonds = [Bond(b, b+1) for b in range(self.nbonds)]
        self.angles = [Angle(a, a+1, a+2) for a in range(self.nangles)]
        self.dihedrals = [Dihedral(d, d+1, d+2, d+3) for d in range(self.ndihedrals)]

        # Call the initialization functions
        self.gen_coords()
        self.gen_types()
        self.gen_bonds()
        self.gen_angles()

    def gen_coords(self):
        """Generate the coordinates of atoms in the molecule

        Returns:
            None
        """
        bondlength = Trappe['bondlength']
        bondangle = Trappe['angles'][AngleType.BACKBONE]['theta0'] * pi / 180.
        # Assign coordinates
        for atom in range(0, self.natoms, 2):
            self.atoms[atom].coords[0] = atom * bondlength * sin(bondangle / 2)
            self.atoms[atom].coords[1] = 0.0
            self.atoms[atom].coords[2] = 0.0
        for atom in range(1, self.natoms, 2):
            self.atoms[atom].coords[0] = atom * bondlength * sin(bondangle / 2)
            self.atoms[atom].coords[1] = cos(bondangle / 2)
            self.atoms[atom].coords[2] = 0.0

        return

    def gen_types(self):
        """
        Assign the types to each atom in the molecule
        """
        self.atoms[0].type = AtomType.CH3
        for atom in range(1, self.natoms):
            self.atoms[atom].type = AtomType.CH2
        self.atoms[-1].type = AtomType.CH3

    def gen_bonds(self):
        """Assign the bonds and bond-type to the molecule

        """
        for bond in range(0, self.nbonds):
            self.bonds[bond].type = BondType.SAT
            self.bonds[bond].atoms = (bond + 0, bond + 1)
        return

    def gen_angles(self):
        """Assign the angles and angle-type to the molecule

        """
        for angle in range(0, self.nangles):
            self.angles[angle].type = AngleType.BACKBONE
            self.angles[angle].atoms = (angle + 0, angle + 1, angle + 2)

    def gen_adihedrals(self):
        """Assign the dihedrals and dihedral-type to the molecule

        """
        for dihedral in range(0, self.ndihedrals):
            self.dihedrals[dihedral].type = DihedralType.BACKBONE
            self.dihedrals[dihedral].atoms = (dihedral + 0, dihedral + 1, dihedral + 2, dihedral + 3)

    def set_indices(self, molidx):
        """Set index of molecule and assign ID to members (atoms, bonds etc.) for LAMMPS config

        Parameters:
            molidx (int): Index of molecule
        """
        self.index = molidx
        for idx in range(self.natoms):
            self.atoms[idx].index = molidx * self.natoms + idx
        for idx in range(self.nbonds):
            self.bonds[idx].index = molidx * self.nbonds + idx
        for idx in range(self.nangles):
            self.angles[idx].index = molidx * self.nangles + idx
        for idx in range(self.ndihedrals):
            self.dihedrals[idx].index = molidx * self.ndihedrals + idx

        return

    def print_coords(self):
        """
        Print the coordinates of atoms in the molecule for debugging
        :return: None
        """
        for atom in range(self.natoms):
            print('Atom:{:d} Type:{:s} {:.3f} {:.3f} {:.3f}'.format(atom, self.atoms[atom].type,
                                                                    self.atoms[atom].coords[0],
                                                                    self.atoms[atom].coords[1],
                                                                    self.atoms[atom].coords[2]))

    def translate(self, delta: (float, float, float)):
        """
        Translate the molecule by given distance in x,y,z directions
        :return: A copy of the translated coordinates
        """
        for atom in range(self.natoms):
            self.atoms[atom].coords += delta
        return

    def move_to_coords(self, pos: (float, float, float)) -> None:
        """
        Move the molecule to the given position
        """
        self.translate(pos - self.base_coords)
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
        nvec = np.cross(self.paxis, axis)  # normal vector to current and target axis
        theta = np.arcsin(np.linalg.norm(nvec))  # angle between current and target axis
        norm = np.linalg.norm(nvec)
        if norm > 0.0:
            nvec = nvec / norm  # normalize
        r = Rotation.from_rotvec(theta * nvec)  # get the rotation matrix
        for atom in range(self.natoms):
            self.atoms[atom].coords = r.apply(self.atoms[atom].coords)  # apply rot to each atom
        self.paxis = axis
        return
