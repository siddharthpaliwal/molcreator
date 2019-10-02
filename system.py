"""
Contains the main system definitions
"""
import os
from molecule import Molecule
from geometry import Planar
from trappe import Trappe


def write_settings(path: str) -> object:
    """
    Write the force-field settings to LAMMPS input file
    :type path: String
    :rtype: object
    :param path: Path to the directory where the file should be placed
    """
    folder = os.path.abspath(path)
    with open(folder + '/system.in.settings', 'w') as f:
        f.write("""
pair_coeff 1 1 lj/charmm/coul/charmm 0.091411522 3.95
pair_coeff 2 2 lj/charmm/coul/charmm 0.194746286 3.75
pair_coeff 3 3 lj/charmm/coul/charmm 0.294106636 3.73
bond_coeff     1    harmonic   120.0   1.54
angle_coeff    1    harmonic   62.0022 114 
dihedral_coeff 1 opls 1.411036 -0.271016 3.145034 0.0 
group TraPPE type 1 2 3 

pair_coeff 4 4 lj/charmm/coul/charmm 0.4610313512 3.62
pair_coeff 5 5 lj/charmm/coul/charmm 0.0 0.0 
bond_coeff     2   harmonic   120.0   1.82
bond_coeff     3    harmonic   120.0   1.34
angle_coeff    2        harmonic   62.0022 114.0
angle_coeff    3        harmonic   33.6135 96.0
dihedral_coeff 2 opls -0.20686795 0.0733675754 1.2175996962 0.0    
                """)
    return


class System(object):
    """
    Define a LAMMPS system object with all the necessary attributes
    """
    tolerance = 2.0

    def __init__(self,
                 nmol: int,
                 box: (float, float, float),
                 moltype,
                 natompermol,
                 origin=(0.0, 0.0, 1.0),
                 normal=(0.0, 0.0, 1.0),
                 geom_type='planar'):
        """
        Create a system with particles allocated on the manifold
        and default attributes for LAMMPS
        Args:
            nmol: No. of molecules
            origin: Box origin
            box: box dimensions
            moltype: string name of alkane
            natompermol: No. of atoms per molecule
            normal: vector normal to the planar geometry
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
        assert (len(moltype) == 3 and isinstance(moltype, str))
        assert (isinstance(natompermol, int))

        self.nmol = nmol
        self.natompermol = natompermol
        self.origin = origin
        self.box = box

        self.geometry = self.gen_manifold(geom_type, normal)
        self.molecules = self.gen_molecules(moltype, natompermol)
        self.atoms = self.nmol * self.natompermol

    def gen_manifold(self, type, normal):
        """
        Generate the manifold and seeds
        :param type: Type of manifold ('planar' or 'sphere')
        :param normal: Direction of normal (in case of 'planar')
        :return:
        """
        geometry = None
        if type == 'planar':
            print("Generating a planar manifold")
            geometry = Planar(self.nmol, boxlen=self.box, normal=normal)
            print("Generating seeds")
            geometry.generate_seeds()
            print("Done")
        return geometry

    def gen_molecules(self, moltype, natompermol) -> list:
        """
        Generate 'nmol' molecules in the system and give them correct coordinates
        :param moltype: String defining the type of linear alkane molecule
        :param natompermol: No. of atoms per molecule
        :return: List: list of generated molecules
        """
        molecules = [Molecule(moltype, natompermol) for _ in range(self.nmol)]
        if self.geometry is not None:
            for i in range(self.nmol):
                molecules[i].rotate_paxis(self.geometry.normal)
                molecules[i].move_to_coords(self.geometry.seeds[i])
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

    def write_coords_lmp(self, path: str) -> object:
        """
        Write the generated atom coordinates to a LAMMPS data file
        :type path: String
        :rtype: object
        :param path: Path to the directory where the file should be placed
        """
        folder = os.path.abspath(path)
        with open(folder + '/system.data', 'w') as f:
            f.write("LAMMPS Description\n")
            f.write("\n")
            f.write("{:d}  atoms\n".format(self.natompermol * self.nmol))
            f.write("{:d}  bonds\n".format(len(self.molecules[0].bond) * self.nmol))
            f.write("{:d}  angles\n".format(0))
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
                f.write("{:d}  {:f}\n".format(itype+1, mass, atomtype))
            f.write("\n")

            # Atom coordinates
            f.write("Atoms\n\n")
            atom_count = 0
            for imol, mol in enumerate(self.molecules):
                for iatom in range(mol.natoms):
                    atom_count += 1
                    f.write("{:d} {:d} {:d} {:.3f} {:.5f} {:.5f} {:.5f}\n".format(atom_count, imol+1, 1, 0.0,
                                                                                  mol.atom_coords[iatom, 0],
                                                                                  mol.atom_coords[iatom, 1],
                                                                                  mol.atom_coords[iatom, 2]))

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
