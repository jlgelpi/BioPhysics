#!/usr/bin/env python3
"""Create a vdwprm-style file and a coordinate file with AMBER types.

Usage: python prepare_amber_files.py [topology] [coordinates]

The coordinate file must be the GRO produced with the expanded topology.
The PDBQT-like output contains AMBER types, not AutoDock types. Charges retain
four decimal places. Sigma is converted from nm to Angstrom and epsilon from
kJ/mol to kcal/mol. Surface coefficients are approximate categories taken
from the reference vdwprm, not parameters supplied by the AMBER force field.
Only type-level Lennard-Jones parameters are exported, not pair overrides or
1-4 interactions.
"""

import argparse
import math
import re
from pathlib import Path

from GROMACS.extract_atom_parameters import read_topology

HERE = Path(__file__).resolve().parent


def read_type_metadata(path):
    metadata = {}
    section = None
    with path.open(encoding="utf-8") as source:
        for line in source:
            line = line.split(";", 1)[0].strip()
            if line.startswith("[") and line.endswith("]"):
                section = line[1:-1].strip().lower()
            elif line and section == "defaults":
                if line.split()[:2] != ["1", "2"]:
                    raise ValueError("only Lennard-Jones combination rule 2 is supported")
            elif line and section == "atomtypes":
                fields = line.split()
                metadata[fields[0]] = (int(fields[1]), float(fields[2]))
    return metadata


def read_surface_coefficients(path):
    coefficients = {}
    with path.open(encoding="utf-8") as source:
        for line in source:
            fields = line.split()
            if not fields or fields[0].startswith("#") or fields[0] == "TIPXEP":
                continue
            coefficients[fields[0]] = float(fields[4])
    return coefficients


def surface_category(atom_type, atomic_number):
    if atomic_number == 1:
        if atom_type in {"H", "HS"}:
            return "HN"
        if atom_type == "HO":
            return "HO"
        if atom_type in {"HW", "HW_spc", "HW_tip4p", "HW_tip4pew"}:
            return "HW"
        if atom_type in {"HC", "H1", "H2", "H3", "H4", "H5", "HA", "HP"}:
            return "H"
    elif atomic_number == 8 and atom_type.startswith("OW"):
        return "OW"
    elif atomic_number in {6, 7, 8, 15, 16}:
        return {6: "C", 7: "N", 8: "O", 15: "P", 16: "S"}[atomic_number]
    raise ValueError(f"no reference surface category for atom type {atom_type}")


def expanded_atoms(molecule_types, molecules):
    for molecule, count in molecules:
        match = re.fullmatch(r".*_chain_([A-Za-z0-9])", molecule)
        chain = match.group(1) if match else " "
        for _ in range(count):
            for atom in molecule_types[molecule]:
                yield chain, atom


def make_outputs(topology_path, gro_path, reference_path):
    atomtypes, molecule_types, molecules = read_topology(topology_path)
    metadata = read_type_metadata(topology_path)
    coefficients = read_surface_coefficients(reference_path)
    atoms = list(expanded_atoms(molecule_types, molecules))
    used_types = sorted({atom[1][1] for atom in atoms})
    parameters = ["# type epsilon(kcal/mol) sigma(A) mass fsrf\n"]
    for name in used_types:
        atomic_number, mass = metadata[name]
        category = surface_category(name, atomic_number)
        if category not in coefficients or not math.isfinite(coefficients[category]):
            raise ValueError(f"missing or invalid surface coefficient for {category}")
        sigma, epsilon = atomtypes[name]
        parameters.append(
            f"{name:<8s} {float(epsilon) / 4.184:.4f} "
            f"{float(sigma) * 10:.4f} {mass:.5f} {coefficients[category]:.6f}\n"
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
    parser.add_argument("topology", nargs="?", type=Path, default=HERE / "2vgb.top.txt")
    parser.add_argument("coordinates", nargs="?", type=Path, default=HERE / "2vbg_fixed.gro")
    parser.add_argument("--reference", type=Path, default=HERE.parent / "data" / "vdwprm")
    parser.add_argument("--parameters", type=Path, default=HERE / "2vgb_amber.vdwprm")
    parser.add_argument("--pdbqt", type=Path, default=HERE / "2vgb_amber.pdbqt")
    args = parser.parse_args()
    try:
        outputs = (args.parameters, args.pdbqt)
        if outputs[0].resolve() == outputs[1].resolve():
            raise ValueError("parameter and PDBQT outputs must be different files")
        for path in outputs:
            if path.exists():
                raise ValueError(f"output already exists: {path}")
        parameters, records = make_outputs(args.topology, args.coordinates, args.reference)
        with args.parameters.open("x", encoding="utf-8") as output:
            output.writelines(parameters)
        with args.pdbqt.open("x", encoding="utf-8") as output:
            output.writelines(records)
    except (OSError, ValueError, KeyError) as exc:
        parser.exit(1, f"error: {exc}\n")
    print(f"Wrote {len(parameters) - 1} types to {args.parameters}")
    print(f"Wrote {len(records)} atoms to {args.pdbqt}")


if __name__ == "__main__":
    main()
