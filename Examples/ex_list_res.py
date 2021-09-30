#!/usr/bin/env python
#
""" Simple program to print ARG residues iteration over atoms """

from Bio.PDB.PDBParser import PDBParser

parser = PDBParser()

st = parser.get_structure('1UBQ', '1ubq.pdb')

selec = []
aa = ["ARG"]

for at in st.get_atoms():
    if at.get_parent().get_resname() in aa:
        selec.append(at)

print("Coordinates:")
for atom in selec:
    print(atom.get_parent().get_resname(), atom.get_parent().id,
          atom.get_name(), atom.get_coord())
