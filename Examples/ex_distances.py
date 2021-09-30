#!/usr/bin/env python
#
""" Simple program to search contacts """

from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser

MAXDIST = 20  # Define distance for a  contact

parser = PDBParser(PERMISSIVE=1)

# load structure from PDB file

st = parser.get_structure('1UBQ', '1ubq.pdb')

select = []

#Select only CA atoms

for at in st.get_atoms():
    if at.id == 'CA':
        select.append(at)
        print("ATOM:", at.get_parent().get_resname(),at.get_parent().id[1], at.id)	

# Preparing search
nbsearch = NeighborSearch(select)

print("NBSEARCH:")

#Searching for contacts under HBLNK

ncontact = 1

for at1, at2 in nbsearch.search_all(MAXDIST):
    print("Contact: ", ncontact)
    print("at1", at1, at1.get_serial_number(), at1.get_parent().get_resname())
    print("at2", at2, at2.get_serial_number(), at2.get_parent().get_resname())
    print()
    ncontact += 1
