"""
Contains the main system definitions
"""
import os
from molcreator.molecule import Molecule
from molcreator.geometry import Planar
from molcreator.trappe import Trappe


class System(object):
    """
    Define a LAMMPS system object with all the necessary attributes
    """
    tolerance = 2.0
    force_field = 'TraPPE-UA'
    units = 'real'
    atom_style = 'full'
    boundary = 'p p p'  # default
    bond_style = 'hybrid harmonic'
    angle_style = 'hybrid harmonic'
    dihedral_style = 'hybrid harmonic'
    pair_style = 'hybrid lj/charmm/coul/charmm'
    pair_modify = 'mix arithmetic'  # mixing rule for pairwise interactions
    special_bonds = 'lj 0.0 0.0 0.0'  # any special 1-3 1-4 interactions

    def __init__(self,
                 molecule: Molecule,
                 nmol: int,
                 box: (float, float, float),
                 origin=(0.0, 0.0, 0.0)):
        """Create a system with particles allocated on the manifold
        and default attributes for LAMMPS

        Args:
            nmol: No. of molecules
            origin: Box origin
            box: box dimensions
            molecule (Molecule): A molecule object to replicate
        """
        assert (nmol < 10000)
        self.nmol = nmol
        self.origin = origin
        self.box = box

        self.geometry = None
        self.molecules = None
        self.natoms = self.nmol * molecule.natoms

    def gen_manifold(self, gtype, normal, tol) -> None:
        """Generate the manifold and seeds

        Args:
            gtype (str): Type of manifold ('planar' or 'sphere')
            normal (float, float, float): Direction of normal (in case of 'planar')
            tol (float): Tolerance distance between molecules

        Returns:
            None
        """
        geometry = None
        if gtype == 'planar':
            print("Generating a planar manifold")
            geometry = Planar(self.nmol, boxlen=self.box, normal=normal)
            print("Generating seeds")
            # geometry.generate_seeds()
            # geometry.generate_seeds_poisson(tol=tol)
            geometry.generate_seeds_square(tol)
            print("Done")
        else:
            print('other geometries not implemented yet!')
        self.geometry = geometry
        return

    def gen_molecules(self, molecule: Molecule) -> None:
        """Generate 'nmol' molecules in the system and give them correct coordinates

        Args:
            molecule: A :class Molecule object to replicate

        Returns:
            None
        """
        self.molecules = [Molecule.from_molecule(molecule, idx) for idx in range(self.nmol)]
        if self.geometry is not None:
            atom_count = 0
            for molidx, molecule in enumerate(self.molecules):
                atom_count += molecule.natoms
                molecule.rotate_paxis(self.geometry.normal)
                molecule.move_to_coords(self.geometry.seeds[molidx])
                molecule.set_indices(molidx)
        else:
            print("Not implemented yet!")
            pass
        return

    @classmethod
    def write_settings(cls, path: str) -> object:
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
            f.write("{:d}  dihedrals\n".format(self.molecules[0].ndihedrals * self.nmol))
            f.write("{:d}  impropers\n".format(0))
            f.write("\n")
            f.write("{:d}  atom types\n".format(len(Trappe['masses'])))
            f.write("{:d}  bond types\n".format(1))
            f.write("{:d}  angle types\n".format(1))
            f.write("{:d}  dihedral types\n".format(1))
            f.write("\n")
            f.write("{:.5f} {:.5f} xlo xhi\n".format(self.origin[0], self.box[0] + self.origin[0]))
            f.write("{:.5f} {:.5f} ylo yhi\n".format(self.origin[1], self.box[1] + self.origin[1]))
            f.write("{:.5f} {:.5f} zlo zhi\n".format(self.origin[2], self.box[2] + self.origin[2]))
            f.write("\n")

            # Masses
            f.write("Masses\n\n")
            for itype, (atomtype, mass) in enumerate(Trappe['masses'].items()):
                f.write("{:d}  {:f}\n".format(itype + 1, mass, atomtype))
            f.write("\n")

            # Atom data
            f.write("Atoms\n\n")
            for mol in self.molecules:
                for atom in mol.atoms:
                    f.write(
                        f'{atom.index + 1:d} {mol.index + 1:d} {atom.type:d} {0.0:.3f} {atom.coords[0]:.5f}'
                        f' {atom.coords[1]:.5f} {atom.coords[2]:.5f}\n')
            f.write("\n")

            # Bond data
            f.write("Bonds\n\n")
            for mol in self.molecules:
                for bond in mol.bonds:
                    f.write(f'{bond.index + 1:d} {bond.type:d} {bond.atoms[0] + 1:d} {bond.atoms[1] + 1:d}\n')
            f.write("\n")

            # Angles data
            f.write("Angles\n\n")
            for mol in self.molecules:
                for angle in mol.angles:
                    f.write(
                        f'{angle.index + 1:d} {angle.type:d} {angle.atoms[0] + 1:d} {angle.atoms[1] + 1:d}'
                        f' {angle.atoms[2] + 1:d}\n')
            f.write("\n")

            # Dihedrals data
            f.write("Dihedrals\n\n")
            for mol in self.molecules:
                for dihedral in mol.dihedrals:
                    f.write(
                        f'{dihedral.index + 1:d} {dihedral.type:d} {dihedral.atoms[0] + 1:d} {dihedral.atoms[1] + 1:d}'
                        f' {dihedral.atoms[2] + 1:d} {dihedral.atoms[3] + 1:d}\n')
            f.write("\n")
