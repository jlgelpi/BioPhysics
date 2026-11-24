#!/usr/bin/env python3
"""Create a vdwprm-style file and a coordinate file with AMBER types from GROMACS Topology.

Usage: python prepare_amber_files.py [topology] [coordinates] [--reference REFERENCE] [--parameters PARAMETERS] [--pdbqt PDBQT]

The coordinate file must be the GRO produced with the expanded topology.
The PDBQT-like output contains Forcefield types. Charges retain
four decimal places. Sigma is converted from nm to Angstrom and epsilon from
kJ/mol to kcal/mol. Surface coefficients are approximate categories taken
from the reference file, not parameters supplied by the AMBER force field.
Only type-level Lennard-Jones parameters are exported, not pair overrides or
1-4 interactions.
"""

import argparse
import csv
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

SURFACE_COEFFICIENTS = {
    "A": 0.111, # Aromatic C
    "C": 0.019, # Aliphatic C
    "N": -0.124, 
    "P": 0.000,
    "S": 0.026,
    "H": 0.0, # All hydrogen types 
    "OH": -0.043, # Hydroxyl O
    "OC": -0.069, # Carboxylate O
    "O": -0.031,  # Carbonyl O
}

SURFACE_MAP_TYPES = {
    "CT": "C",
    "C*": "A",
    "CB": "A",
    "CN": "A",
    "CW": "A",
    "CA": "A",
    "CR": "A",
    "CC": "A",
    "CV": "A",
    "O2": "OC",
    "OH": "OH",
    "SH": "S",
    "NA": "N",
    "NB": "N",
    "N2": "N",
    "N3": "N",
    "HC": "H",    
    "H1": "H",
    "H4": "H",
    "H5": "H",
    "HA": "H",
    "HO": "H",
    "HP": "H",
    "HS": "H",
}    

def read_topology(path):
    atomtypes = {}
    type_metadata = {}
    molecule_types = {}
    molecules = []
    section = None
    current_molecule = None
    defaults = None

    with path.open(encoding="utf-8") as topology:
        for line_number, line in enumerate(topology, 1):
            line = line.split(";", 1)[0].strip()
            if not line:
                continue
            if line.startswith("#"):
                raise ValueError(f"{path}:{line_number}: preprocessing directives are not supported")
            if line.startswith("[") and line.endswith("]"):
                section = line[1:-1].strip().lower()
                continue
            fields = line.split()
            try:
                if section == "defaults":
                    if defaults is not None or len(fields) < 2:
                        raise ValueError("expected one [ defaults ] entry")
                    defaults = (int(fields[0]), int(fields[1]))
                    if defaults[0] != 1 or defaults[1] not in (2, 3):
                        raise ValueError("only Lennard-Jones combination rules 2 and 3 are supported")
                elif section == "atomtypes":
                    if len(fields) != 7:
                        raise ValueError("expected name, atomic number, mass, charge, ptype, sigma, epsilon")
                    name = fields[0]
                    if name in atomtypes:
                        raise ValueError(f"duplicate atom type {name}")
                    atomic_number, mass = int(fields[1]), float(fields[2])
                    sigma, epsilon = fields[-2:]
                    if not all(math.isfinite(float(value)) and float(value) >= 0
                               for value in (sigma, epsilon)):
                        raise ValueError(f"invalid Lennard-Jones parameters for {name}")
                    atomtypes[name] = (sigma, epsilon)
                    type_metadata[name] = (atomic_number, mass)
                elif section == "moleculetype":
                    if len(fields) < 2:
                        raise ValueError("invalid molecule type")
                    current_molecule = fields[0]
                    if current_molecule in molecule_types:
                        raise ValueError(f"duplicate molecule type {current_molecule}")
                    molecule_types[current_molecule] = []
                elif section == "atoms":
                    if current_molecule is None or len(fields) < 7:
                        raise ValueError("expected a molecule type and an atom with an explicit charge")
                    nr, atom_type, resnr, residue, atom, _, charge = fields[:7]
                    if int(nr) != len(molecule_types[current_molecule]) + 1:
                        raise ValueError("atom numbers must be consecutive and start at 1")
                    int(resnr)
                    if not math.isfinite(float(charge)):
                        raise ValueError("atom charge must be finite")
                    molecule_types[current_molecule].append(
                        (nr, atom_type, resnr, residue, atom, charge)
                    )
                elif section == "molecules":
                    if len(fields) != 2 or int(fields[1]) < 0:
                        raise ValueError("expected a molecule name and non-negative count")
                    molecules.append((fields[0], int(fields[1])))
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number}: {exc}") from exc

    if defaults is None or not atomtypes or not molecules:
        raise ValueError(f"{path}: missing defaults, atom types, or system molecules")
    for name, count in molecules:
        if name not in molecule_types or not molecule_types[name]:
            raise ValueError(f"{path}: missing atoms for molecule type {name}")
        if count:
            for _, atom_type, _, _, _, _ in molecule_types[name]:
                if atom_type not in atomtypes:
                    raise ValueError(f"{path}: undefined atom type {atom_type} in {name}")
    return atomtypes, type_metadata, molecule_types, molecules

def expanded_atoms(molecule_types, molecules):
    for molecule, count in molecules:
        match = re.fullmatch(r".*_chain_([A-Za-z0-9])", molecule)
        chain = match.group(1) if match else " "
        for _ in range(count):
            for atom in molecule_types[molecule]:
                yield chain, atom



def make_outputs(topology_path, gro_path):
    atomtypes, metadata, molecule_types, molecules = read_topology(topology_path)
    coefficients = SURFACE_COEFFICIENTS
    atoms = list(expanded_atoms(molecule_types, molecules))
    used_types = sorted({atom[1][1] for atom in atoms})
    parameters = ["# type epsilon(kcal/mol) sigma(A) mass fsrf\n"]
    for name in used_types:
        atomic_number, mass = metadata[name]
        if name not in SURFACE_COEFFICIENTS:
            coefficients[name] = SURFACE_COEFFICIENTS[SURFACE_MAP_TYPES.get(name, name)]
        else:
            coefficients[name] = SURFACE_COEFFICIENTS[name]
        
        sigma, epsilon = atomtypes[name]
        parameters.append(
            f"{name:<8s} {float(epsilon) / 4.184:.4f} "
            f"{float(sigma) * 10:.4f} {mass:6.3f} {coefficients[name]:6.3f}\n"
        )

    records = []
    with gro_path.open(encoding="utf-8") as gro:
        if not gro.readline():
            raise ValueError(f"{gro_path}: missing GRO title")
        count = int(gro.readline().strip())
        if count != len(atoms):
            raise ValueError(f"{gro_path}: {count} coordinates but {len(atoms)} topology atoms")
        if count > 99999:
            raise ValueError("PDBQT atom serials cannot exceed 99999")
        for index, (chain, atom) in enumerate(atoms, 1):
            _, atom_type, resnr, residue, name, charge = atom
            line = gro.readline()
            if len(line) < 44:
                raise ValueError(f"{gro_path}: incomplete atom record {index}")
            gro_resnr = int(line[:5])
            gro_residue, gro_name = line[5:10].strip(), line[10:15].strip()
            if (gro_resnr, gro_residue, gro_name) != (
                int(resnr) % 100000, residue[:5], name[:5]
            ):
                raise ValueError(
                    f"{gro_path}: atom {index} does not match topology: "
                    f"{gro_resnr} {gro_residue} {gro_name} != {resnr} {residue} {name}"
                )
            if len(residue) > 3 or len(name) > 4 or not 0 <= int(resnr) <= 9999:
                raise ValueError(f"atom {index} cannot be represented in PDBQT fields")
            coords = [float(line[start:start + 8]) * 10 for start in (20, 28, 36)]
            charge_value = float(charge)
            if not all(math.isfinite(value) for value in (*coords, charge_value)):
                raise ValueError(f"atom {index} has non-finite coordinates or charge")
            xyz = "".join(f"{value:8.3f}" for value in coords)
            charge_field = f"{charge_value:7.4f}"
            if len(xyz) != 24 or len(charge_field) != 7:
                raise ValueError(f"atom {index} coordinates or charge exceed PDBQT field width")
            records.append(
                f"ATOM  {index:5d} {name:>4s} {residue:3s} {chain}{int(resnr):4d}    "
                f"{xyz}{1.0:6.2f}{0.0:6.2f}   {charge_field} {atom_type}\n"
            )
        if not gro.readline().strip():
            raise ValueError(f"{gro_path}: missing GRO box")
    return parameters, records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("topology", nargs="?", type=Path)
    parser.add_argument("coordinates", nargs="?", type=Path)
    parser.add_argument("--output_prefix", type=Path)
    parser.add_argument("--not_overwrite", action="store_true")
    args = parser.parse_args()
    args.parameters = args.output_prefix.with_suffix(".vdwprm")
    args.pdbqt = args.output_prefix.with_suffix(".pdbqt")
    try:
        outputs = [args.parameters, args.pdbqt]
        if len({path.resolve() for path in outputs}) != len(outputs):
            raise ValueError("output files must be different")
        for path in outputs:
            if path.exists():
                if args.not_overwrite:
                    raise ValueError(f"Output already exists: {path}")
                print(f"Output already exists: {path}, overwriting...")
        parameters, records = make_outputs(args.topology, args.coordinates)
        with args.parameters.open("w", encoding="utf-8") as output:
            output.writelines(parameters)
        with args.pdbqt.open("w", encoding="utf-8") as output:
            output.writelines(records)
    except (OSError, ValueError, KeyError) as exc:
        parser.exit(1, f"error: {exc}\n")
    print(f"Wrote {len(parameters) - 1} types to {args.parameters}")
    print(f"Wrote {len(records)} atoms to {args.pdbqt}")


if __name__ == "__main__":
    main()
