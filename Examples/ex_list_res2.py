#!/usr/bin/env python
#
""" Simple program to print residue atoms """


from Bio.PDB.PDBParser import PDBParser


parser = PDBParser()

st = parser.get_structure('1UBQ', '1ubq.pdb')

selec = []
res_num = 20

res = st[0]["A"][res_num]

print ("Residue", res.get_resname(), res.id[1])
print("Atoms:")
for atom in res.get_atoms(): 
    print(res.get_resname(), res.id,
              atom.get_name(), atom.get_coord())
