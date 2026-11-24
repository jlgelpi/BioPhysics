#!/usr/bin/env python3
"""Export atom charges and type-level Lennard-Jones parameters as CSV.

Usage:
    python extract_atom_parameters.py [2vgb_top.txt] > atom_parameters.csv

The input must be an expanded topology with the atomtypes layout used by
2vgb_top.txt and Lennard-Jones combination rule 2 or 3. Sigma is in nm,
epsilon in kJ/mol, and charge in elementary charge units. These are atom-type
parameters, not pair-specific overrides or 1-4 interaction parameters.
For alchemical atoms, the A-state type and charge are exported.
"""

import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path


def read_topology(path):
    atomtypes = {}
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
                    sigma, epsilon = fields[-2:]
                    if not all(math.isfinite(float(value)) and float(value) >= 0
                               for value in (sigma, epsilon)):
                        raise ValueError(f"invalid Lennard-Jones parameters for {name}")
                    atomtypes[name] = (sigma, epsilon)
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
    return atomtypes, molecule_types, molecules


def write_csv(output, atomtypes, molecule_types, molecules):
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow([
        "atom_index", "molecule_type", "molecule_copy", "atom_number",
        "residue_number", "residue_name", "atom_name", "atom_type",
        "charge_e", "sigma_nm", "epsilon_kj_mol",
    ])
    atom_index = 0
    copies = defaultdict(int)
    for molecule_name, count in molecules:
        for _ in range(count):
            copies[molecule_name] += 1
            for nr, atom_type, resnr, residue, atom, charge in molecule_types[molecule_name]:
                atom_index += 1
                sigma, epsilon = atomtypes[atom_type]
                writer.writerow([
                    atom_index, molecule_name, copies[molecule_name], nr,
                    resnr, residue, atom, atom_type, charge, sigma, epsilon,
                ])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "topology", nargs="?", type=Path,
        default=Path(__file__).with_name("2vgb_top.txt"),
        help="expanded topology file (default: 2vgb_top.txt next to this script)",
    )
    args = parser.parse_args()
    try:
        topology = read_topology(args.topology)
        write_csv(sys.stdout, *topology)
    except (OSError, ValueError) as exc:
        parser.exit(1, f"error: {exc}\n")


if __name__ == "__main__":
    main()
