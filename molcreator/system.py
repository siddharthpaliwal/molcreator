"""
Contains the main system definitions
"""
import os
import copy
from typing import List
from molcreator.molecule import Molecule
from molcreator.geometry import Planar
from molcreator.trappe import Trappe


def write_settings(path: str) -> object:
    """Write the force-field settings to LAMMPS input file

    Args:
        path: Path to the directory where the file should be placed

    Returns:
        None
    """
    folder = os.path.abspath(path)
    with open(folder + '/system.in.settings', 'w') as f:
        f.write("")
        f.write("pair_coeff 1 1 lj/charmm/coul/charmm 0.091411522 3.95\n")
        f.write("pair_coeff 2 2 lj/charmm/coul/charmm 0.194746286 3.75\n")
        f.write("pair_coeff 3 3 lj/charmm/coul/charmm 0.294106636 3.73\n")
        f.write("bond_coeff     1    harmonic   120.0   1.54\n")
        f.write("angle_coeff    1    harmonic   62.0022 114\n")
        f.write("dihedral_coeff 1 opls 1.411036 -0.271016 3.145034 0.0 \n")
        f.write("group TraPPE type 1 2 3 \n")
        f.write("pair_coeff 4 4 lj/charmm/coul/charmm 0.4610313512 3.62\n")
        f.write("pair_coeff 5 5 lj/charmm/coul/charmm 0.0 0.0 \n")
        f.write("bond_coeff     2   harmonic   120.0   1.82\n")
        f.write("bond_coeff     3    harmonic   120.0   1.34\n")
        f.write("angle_coeff    2        harmonic   62.0022 114.0\n")
        f.write("angle_coeff    3        harmonic   33.6135 96.0\n")
        f.write("dihedral_coeff 2 opls -0.20686795 0.0733675754 1.2175996962 0.0\n")
    return


class System(object):
    """
    Define a LAMMPS system object with all the necessary attributes
    """
    tolerance = 2.0

    def __init__(self,
                 molecule: Molecule,
                 nmol: int,
                 box: (float, float, float),
                 origin=(0.0, 0.0, 1.0),
                 normal=(0.0, 0.0, 1.0),
                 geom_type='planar'):
        """Create a system with particles allocated on the manifold
        and default attributes for LAMMPS

        Args:
            nmol: No. of molecules
            origin: Box origin
            box: box dimensions
            normal: vector normal to the planar geometry
            molecule (Molecule): A molecule object to replicate in the system
        """
        self.force_field = 'TraPPE-UA'
        self.units = 'real'
        self.atom_style = 'full'
        self.boundary = 'p p p'  # default
        self.bond_style = 'hybrid harmonic'
        self.angle_style = 'hybrid harmonic'
        self.dihedral_style = 'hybrid harmonic'
        self.pair_style = 'hybrid lj/charmm/coul/charmm'
        self.pair_modify = 'mix arithmetic'  # mixing rule for pairwise interactions
        self.special_bonds = 'lj 0.0 0.0 0.0'  # any special 1-3 1-4 interactions

        assert (nmol < 10000)

        self.nmol = nmol
        self.origin = origin
        self.box = box

        self.geometry = self.gen_manifold(geom_type, normal)
        self.molecules = self.gen_molecules(molecule)
        self.natoms = self.nmol * molecule.natoms

    def gen_manifold(self, gtype, normal):
        """Generate the manifold and seeds

        Args:
            gtype: Type of manifold ('planar' or 'sphere')
            normal: Direction of normal (in case of 'planar')

        Returns:
            None
        """
        geometry = None
        if gtype == 'planar':
            print("Generating a planar manifold")
            geometry = Planar(self.nmol, boxlen=self.box, normal=normal)
            print("Generating seeds")
            geometry.generate_seeds()
            print("Done")
        return geometry

    def gen_molecules(self, molecule: Molecule) -> List[Molecule]:
        """Generate 'nmol' molecules in the system and give them correct coordinates

        Args:
            molecule: A :class Molecule object to replicate

        Returns:
            None
        """
        molecules = [copy.deepcopy(molecule) for _ in range(self.nmol)]

        if self.geometry is not None:
            for idx in range(self.nmol):
                molecules[idx].rotate_paxis(self.geometry.normal)
                molecules[idx].move_to_coords(self.geometry.seeds[idx])
                molecules[idx].set_indices(idx)
        else:
            print("Not implemented yet!")
            pass
            # for i in range(self.nmol):
            #     test = True
            #     while (test):
            #         rtest = molecules[i].coords
            #         rtest = np.random.rand(3) * self.box
            #         for j in range(i):
            #             test = Molecule.check_overlap(rtest, self.coords[j])
            #     molecules[i].set_base_coords(rtest)
        return molecules

    def write_coords_lmp(self, path: str) -> None:
        """Write the generated atom coordinates to a LAMMPS data file

        Args:
            path: Path to the directory where the file should be placed

        Returns:
            None
        """
        folder = os.path.abspath(path)
        with open(folder + '/system.data', 'w') as f:
            f.write("LAMMPS Description\n")
            f.write("\n")
            f.write("{:d}  atoms\n".format(self.natoms))
            f.write("{:d}  bonds\n".format(self.molecules[0].nbonds * self.nmol))
            f.write("{:d}  angles\n".format(self.molecules[0].nangles * self.nmol))
            f.write("{:d}  dihedrals\n".format(0))
            f.write("{:d}  impropers\n".format(0))
            f.write("\n")
            f.write("{:d}  atom types\n".format(len(Trappe['masses'])))
            f.write("{:d}  bond types\n".format(1))
            f.write("{:d}  angle types\n".format(0))
            f.write("{:d}  dihedral types\n".format(0))
            f.write("\n")
            f.write("{:.5f} {:.5f} xlo xhi\n".format(self.origin[0], self.box[0] + self.origin[0]))
            f.write("{:.5f} {:.5f} ylo yhi\n".format(self.origin[1], self.box[1] + self.origin[1]))
            f.write("{:.5f} {:.5f} zlo zhi\n".format(self.origin[2], self.box[2] + self.origin[2]))
            f.write("\n")

            # Masses
            f.write("Masses\n\n")
            for itype, (atomtype, mass) in enumerate(Trappe['masses'].items()):
                # f.write("{:d}  {:f} # {:s}\n".format(itype, mass, atomtype))
                f.write("{:d}  {:f}\n".format(itype + 1, mass, atomtype))
            f.write("\n")

            # Atom data
            f.write("Atoms\n\n")
            for mol in self.molecules:
                for atom in mol.atoms:
                    f.write("{:d} {:d} {:d} {:.3f} {:.5f} {:.5f} {:.5f}\n".format(atom.index,
                                                                                  mol.index,
                                                                                  atom.type,
                                                                                  0.0,
                                                                                  atom.coords[0],
                                                                                  atom.coords[1],
                                                                                  atom.coords[2]))
            f.write("\n")

            # Bond data
            f.write("Bonds\n\n")
            for mol in self.molecules:
                for bond in mol.bonds:
                    f.write("{:d} {:d} {:d} {:d}\n".format(bond.index, bond.type, bond.atoms[0], bond.atoms[1]))
            f.write("\n")

            # Angles data
            f.write("Angles\n\n")
            for mol in self.molecules:
                for angle in mol.angles:
                    f.write("{:d} {:d} {:d} {:d}\n".format(angle.index, angle.type, angle.atoms[0], angle.atoms[1],
                                                           angle.atoms[2]))
            f.write("\n")

            # Dihedrals data
            f.write("Dihedrals\n\n")
            for mol in self.molecules:
                for dihedral in mol.dihedrals:
                    f.write("{:d} {:d} {:d} {:d}\n".format(dihedral.index, dihedral.type, dihedral.atoms[0],
                                                           dihedral.atoms[1],
                                                           dihedral.atoms[2], dihedral.atoms[3]))
            f.write("\n")

            # f.write('ITEM: TIMESTEP\n')
            # f.write('{:d}\n'.format(frame.configuration.step))
            # f.write('ITEM: NUMBER OF ATOMS\n')
            # f.write('{:d}\n'.format(frame.particles.N))
            # f.write('ITEM: BOX BOUNDS\n')
            # f.write('{:.5f} {:.5f}\n'.format(-0.5 * frame.configuration.box[0], 0.5 * frame.configuration.box[0]))
            # f.write('{:.5f} {:.5f}\n'.format(-0.5 * frame.configuration.box[1], 0.5 * frame.configuration.box[1]))
            # f.write('{:.5f} {:.5f}\n'.format(-0.5 * frame.configuration.box[2], 0.5 * frame.configuration.box[2]))
            #
            # f.write('ITEM: ATOMS ID X Y Z\n')
            # for p in range(frame.particles.N):
            #     f.write('{:5d}  {:f} {:f} {:f}\n'.format(p, curdata[p, 0], curdata[p, 1], curdata[p, 2]))
