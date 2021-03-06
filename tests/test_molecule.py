import unittest
from molcreator.molecule import Molecule


class TestMolecule(unittest.TestCase):
    def test_molecule(self):
        ethane = Molecule('ETH', 2, 1, 0, 0)
        # ethane.print_coords()
        self.assertEqual(2, ethane.natoms)

        self.assertAlmostEqual(1.292, ethane.atoms[1].coords[0], 2, "x-coord not a multiple of projected bond-length")

        eicosane = Molecule('EIC', 20, 19, 18, 17)
        self.assertEqual(20, eicosane.natoms)

    def test_butane_natoms(self):
        butane = Molecule('BUT', 4, 3, 2, 1)
        self.assertEqual(4, butane.natoms, "Natoms not correct")

    def test_butane_coords(self):
        butane = Molecule('BUT', 4, 3, 2, 1)
        self.assertAlmostEqual(0.0, butane.atoms[2].coords[2], 2, "z-coord not 0")

    def test_rotation(self):
        butane = Molecule('BUT', 4, 3, 2, 14)
        newaxis = [0., 0., 1.]
        butane.print_coords()
        butane.rotate_paxis(newaxis)
        butane.print_coords()
        self.assertListEqual(butane.paxis, newaxis, "Princ. axis not rotated correctly")

    def test_molecule_move(self):
        mol = Molecule('BUT', 4, 3, 2, 1)
        new_pos = [0.5, 1.5, 0.1]
        mol.move_to_coords(new_pos)
        self.assertListEqual(mol.base_coords, new_pos, "Molecule move incorrect")


if __name__ == '__main__':
    unittest.main(verbosity=2)
