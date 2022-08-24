# SIDECHAIN_TORSIONS: contains for all sidechain torsions the
# PDB atoms that define them
_SIDECHAIN_TORSIONS_AMBER = {
    'ALA':
        [],
    'ARG':
        [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"],
         ["CB", "CG", "CD", "NE"], ["CG", "CD", "NE", "CZ"]],
    'ASN':
        [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "OD1"]],
    'ASP':
        [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "OD1"]],
    'CYS':
        [["N", "CA", "CB", "SG"]],
    'GLN':
        [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"],
         ["CB", "CG", "CD", "OE1"]],
    'GLU':
        [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"],
         ["CB", "CG", "CD", "OE1"]],
    'GLY':
        [],
    'HID':
        [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "ND1"]],
    'HIE':
        [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "ND1"]],
    'HIP':
        [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "ND1"]],
    'ILE':
        [["N", "CA", "CB", "CG1"], ["CA", "CB", "CG1", "CD1"]],
    'LEU':
        [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]],
    'LYS':
        [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"],
         ["CB", "CG", "CD", "CE"], ["CG", "CD", "CE", "NZ"]],
    'MET':
        [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "SD"],
         ["CB", "CG", "SD", "CE"]],
    'PHE':
        [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]],
    'PRO':
        [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"], ["CB", "CG", "CD", "N"]],
    'SER':
        [["N", "CA", "CB", "OG"]],
    'THR':
        [["N", "CA", "CB", "OG1"]],
    'TRP':
        [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]],
    'TYR':
        [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]],
    'VAL':
        [["N", "CA", "CB", "CG1"]]}

_SIDECHAIN_TORSIONS_CHARMM = {
    'ALA':
        [],
    'ARG':
        [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"],
         ["CB", "CG", "CD", "NE"], ["CG", "CD", "NE", "CZ"]],
    'ASN':
        [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "OD1"]],
    'ASP':
        [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "OD1"]],
    'CYS':
        [["N", "CA", "CB", "SG"]],
    'GLN':
        [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"],
         ["CB", "CG", "CD", "OE1"]],
    'GLU':
        [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"],
         ["CB", "CG", "CD", "OE1"]],
    'GLY':
        [],
    'HSD':
        [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "ND1"]],
    'HSE':
        [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "ND1"]],
    'HSP':
        [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "ND1"]],
    'ILE':
        [["N", "CA", "CB", "CG1"], ["CA", "CB", "CG1", "CD"]],
    'LEU':
        [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]],
    'LYS':
        [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"],
         ["CB", "CG", "CD", "CE"], ["CG", "CD", "CE", "NZ"]],
    'MET':
        [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "SD"],
         ["CB", "CG", "SD", "CE"]],
    'PHE':
        [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]],
    'PRO':
        [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD"], ["CB", "CG", "CD", "N"]],
    'SER':
        [["N", "CA", "CB", "OG"]],
    'THR':
        [["N", "CA", "CB", "OG1"]],
    'TRP':
        [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]],
    'TYR':
        [["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]],
    'VAL':
        [["N", "CA", "CB", "CG1"]]}

_BB_ATOMS = ('N', 'CA', 'C', 'O', 'H', 'HN', 'HA', 'HA2', '1H', '2H', '3H', 'OXT')

# VMD selection for all backbone atoms
backbone_atoms = 'name ' + ' or name '.join(_BB_ATOMS)

