"""
Definitions of TraPPE force-field from Siepmann group

# This file stores complete LAMMPS data for the TraPPE model of saturated
# hydrocarbon chains.  In this "united-atom" model, each methyl group is
# represented by a single atom.  Forces between "atoms" are taken from the
# TraPPE force-field. (J Phys Chem B, 1998, volume 102, pp.2569-2577)
"""
Trappe = {
    'masses': {
        'CH4': 16.3307, # gm/mol
        'CH3': 15.2507, # gm/mol
        'CH2': 14.1707 # gm/mol
    },
    'angles': {
        'backbone': {
            'ktheta': 62.0022, # kcal/rad^2
            'theta0': 114 # deg
        }
    },
    'dihedrals': {
        'backbone': {
            'type': 'opls', # https://lammps.sandia.gov/doc/dihedral_opls.html
            'c1': 1.411036, # kcal/mol
            'c2': -0.271016, # kcal/mol
            'c3': 3.145034, # kcal/mol
            'c4': 0.0  # kcal/mol
        }
    },
    'rcut': [12.0,14.0], # inner and outer cutoff for pair-style: lj/charmm/coul/charmm
    'bondlength': 1.540, # Angstroms
}
