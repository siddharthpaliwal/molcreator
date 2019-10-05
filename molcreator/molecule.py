"""
The :class:`~.Molecule` defines the basic unit of a molecule in the system.
The attributes are defined to correspond to a United Atom force field model
for a linear chain alkane.
"""

import numpy as np
from math import sin, cos, pi
import copy
from scipy.spatial.transform import Rotation
from molcreator.trappe import Trappe, AtomType, BondType, AngleType, DihedralType


class Atom(object):
    """Stores the information of each atom

    Attributes:
        type (:class:`molcreator.trappe.AtomType`): Type of Atom
        coords (float, float, float): coordinates of the atom
        index (int): index
    """
    _atom_count = 0

    def __init__(self, coords=None, atype=AtomType.CH2, index: int = 0):
        if coords is None:
            coords = [0.0, 0.0, 0.0]
        Atom._atom_count += 1
        self.type = AtomType(atype)
        self.coords = coords
        self.index = index


class Bond(object):
    """Stores the information of each bond

    Attributes:
        type (:class:`molcreator.trappe.BondType`): Type of Bond
        atoms (int, int): tuple of 2 consecutive atoms connected by the bond
        index (int): index
    """

    def __init__(self, atomi: int, atomj: int, btype: BondType = BondType.SAT):
        """

        Args:
            btype (object): BondType Enum
        """
        self.type = BondType(btype)
        self.atoms = [atomi, atomj]
        self.index = 0


class Angle(object):
    """Stores the information of each angle

    Attributes:
        type (:class:`molcreator.trappe.AngleType`): Type of Angle
        atoms (int, int, int): tuple of 3 consecutive atoms forming the angle
        index (int): index
    """

    def __init__(self, atomi: int, atomj: int, atomk: int, atype=AngleType.BACKBONE):
        self.type = AngleType(atype)
        self.atoms = [atomi, atomj, atomk]
        self.index = 0


class Dihedral(object):
    """Stores the information of each dihedral

    Attributes:
        type (:class:`molcreator.trappe.DihedralType`): Type of Dihedral angle
        atoms (int, int, int, int): tuple of 4 consecutive atoms forming the dihedral angle
        index (int): index
    """

    def __init__(self, atomi: int, atomj: int, atomk: int, atoml: int, atype=DihedralType.BACKBONE):
        self.type = DihedralType(atype)
        self.atoms = [atomi, atomj, atomk, atoml]
        self.index = 0


class Molecule(object):
    R"""Defines the molecule class

    Used to instantiate molecule objects
    that could be modified in space and arranged in the desired manner
    specified by a geometrical manifold

    Args:
        name (string): Chemical/Ref. name of the molecule
        natoms (integer): No. of beads in the molecule representation
        nbonds (integer): No. of beads in the molecule representation
        nangles (integer): No. of beads in the molecule representation
        ndihedrals (integer): No. of beads in the molecule representation
        idx (integer): Index number of the molecule (optional)

    Attributes:
        base_coords: coordinates of the base
        paxis: principle axis unit vector in 3D
        index: index number of the molecule
        atoms (list): list of :class:`~.Atom` objects
        bonds (list): list of :class:`~.Bond` objects
        angles (list): list of :class:`~.Angle` objects
        dihedrals (list): list of :class:`~.Dihedral` objects
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
        self.atoms = [Atom() for _ in range(self.natoms)]
        self.bonds = [Bond(b, b + 1) for b in range(self.nbonds)]
        self.angles = [Angle(a, a + 1, a + 2) for a in range(self.nangles)]
        self.dihedrals = [Dihedral(d, d + 1, d + 2, d + 3) for d in range(self.ndihedrals)]

        # Call the initialization functions
        self.gen_types()
        self.gen_coords()
        self.gen_bonds()
        self.gen_angles()

    @classmethod
    def from_molecule(cls, template_mol, index):
        """Copy constructor

        Args:
            template_mol: Template molecule to copy
            index: current molecule index

        Returns:

        """
        mol = copy.deepcopy(template_mol)
        mol.set_indices(index)
        return mol

    def gen_coords(self):
        """Generate the coordinates of atoms in the molecule

        Returns:
            None
        """
        bondlength = Trappe['bondlength']
        bondangle = Trappe['angles'][AngleType.BACKBONE]['theta0'] * pi / 180.
        # Assign coordinates
        for i, atom in enumerate(self.atoms):
            atom.coords[0] = i * bondlength * sin(bondangle / 2)
            atom.coords[1] = 0.0 if i % 2 == 0 else cos(bondangle / 2)
            atom.coords[2] = 0.0
        return

    def gen_types(self):
        """
        Assign the types to each atom in the molecule
        """
        self.atoms[0].type = AtomType.CH3
        for atom in self.atoms[1:-1]:
            atom.type = AtomType.CH2
        self.atoms[-1].type = AtomType.CH3

    def gen_bonds(self):
        """Assign the bonds and bond-type to the molecule

        """
        for idx, bond in enumerate(self.bonds):
            bond.type = BondType.SAT
            bond.atoms = (idx + 0, idx + 1)
        return

    def gen_angles(self):
        """Assign the angles and angle-type to the molecule

        """
        for idx, angle in enumerate(self.angles):
            angle.type = AngleType.BACKBONE
            angle.atoms = (idx + 0, idx + 1, idx + 2)

    def gen_adihedrals(self):
        """Assign the dihedrals and dihedral-type to the molecule

        """
        for idx, dihedral in enumerate(self.dihedrals):
            dihedral.type = DihedralType.BACKBONE
            dihedral.atoms = (idx + 0, idx + 1, idx + 2, idx + 3)

    def set_indices(self, molidx):
        """Set index of molecule and assign ID to members (atoms, bonds etc.) for LAMMPS config

        Parameters:
            molidx (int): Index of molecule
        """
        self.index = molidx
        for idx, atom in enumerate(self.atoms):
            atom.index = molidx * self.natoms + idx
        for idx, bond in enumerate(self.bonds):
            bond.index = molidx * self.nbonds + idx
            bond.atoms = [self.atoms[idx].index, self.atoms[idx+1].index]
        for idx, angle in enumerate(self.angles):
            angle.index = molidx * self.nangles + idx
            angle.atoms = [self.atoms[idx].index, self.atoms[idx+1].index, self.atoms[idx+2].index]
        for idx, dihedral in enumerate(self.dihedrals):
            dihedral.index = molidx * self.ndihedrals + idx
            dihedral.atoms = [self.atoms[idx].index, self.atoms[idx+1].index, self.atoms[idx+2].index,
                              self.atoms[idx+3].index]
        return

    def print_coords(self):
        """
        Print the coordinates of atoms in the molecule for debugging
        :return: None
        """
        for atom in self.atoms:
            print('Atom:{:d} Type:{:d} {:.3f} {:.3f} {:.3f}'.format(atom, atom.type,
                                                                    atom.coords[0],
                                                                    atom.coords[1],
                                                                    atom.coords[2]))

    def translate(self, delta: (float, float, float)):
        """Translate the molecule by given distance in x,y,z directions

        Args:
            delta (float, float, float): displacement in 3D

        Returns:
            None
        """
        for atom in self.atoms:
            atom.coords += delta
        return

    def move_to_coords(self, pos: (float, float, float)) -> None:
        """Move the molecule to the given position

        Args:
            pos (float, float, float): Target position of the base of the molecule
        """
        self.translate(np.asarray(pos) - np.asarray(self.base_coords))
        self.base_coords = pos
        return

    def get_distance(self, moltest):
        """Calculate the min. distance between any two atoms of self from moltest
        """
        return

    def check_overlap(self, ri, rj):
        """Check if the two molecules overlap
        """
        return True

    def rotate_paxis(self, axis) -> None:
        """rotate the molecule such that the principle axis is now pointing in 'axis' direction

        Args:
            axis (float, float, float): Vector where the :class:`~.Molecule.paxis` should point to
        """
        nvec = np.cross(self.paxis, axis)  # normal vector to current and target axis
        theta = np.arcsin(np.linalg.norm(nvec))  # angle between current and target axis
        norm = np.linalg.norm(nvec)
        if norm > 0.0:
            nvec = nvec / norm  # normalize
        r = Rotation.from_rotvec(theta * nvec)  # get the rotation matrix
        for atom in self.atoms:
            atom.coords = r.apply(atom.coords)  # apply rot to each atom
        self.paxis = axis
        return
