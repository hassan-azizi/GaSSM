#!/usr/bin/env python3
"""Generate a GaSSM input file from a CIF structure and force-field JSON."""

from __future__ import annotations

import argparse
import json
import re
import shlex
import sys
import tomllib
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Atom:
    name: str
    frac_x: float
    frac_y: float
    frac_z: float
    charge: float


@dataclass
class CifData:
    cell_a: float
    cell_b: float
    cell_c: float
    alpha: float
    beta: float
    gamma: float
    atoms: list[Atom]


@dataclass
class AdsorbateSite:
    name: str
    x: float
    y: float
    z: float


class InputError(Exception):
    pass


CELL_KEYS = {
    "_cell_length_a": "cell_a",
    "_cell_length_b": "cell_b",
    "_cell_length_c": "cell_c",
    "_cell_angle_alpha": "alpha",
    "_cell_angle_beta": "beta",
    "_cell_angle_gamma": "gamma",
}


def parse_number(value: str) -> float:
    value = value.strip().strip("'\"")
    if value in {"?", "."}:
        raise ValueError(f"missing numeric value {value!r}")
    value = re.sub(r"\([^)]*\)$", "", value)
    return float(value)


def resolve_path(path_value: str, config_dir: Path) -> Path:
    path = Path(path_value)
    if path.is_absolute():
        return path
    return config_dir / path


def resolve_structure(path_value: str, config_dir: Path) -> Path:
    candidates: list[Path] = []
    raw = Path(path_value)

    if raw.is_absolute():
        candidates.append(raw)
        if raw.suffix != ".cif":
            candidates.append(raw.with_suffix(".cif"))
    else:
        for prefix in (config_dir, config_dir / "Structures"):
            candidates.append(prefix / raw)
            if raw.suffix != ".cif":
                candidates.append(prefix / raw.with_suffix(".cif"))

    for candidate in candidates:
        if candidate.exists():
            return candidate

    tried = "\n  ".join(str(candidate) for candidate in candidates)
    raise InputError(f"Could not find CIF structure {path_value!r}. Tried:\n  {tried}")


def resolve_adsorbate_species(species_value: str, config_dir: Path) -> Path:
    candidates: list[Path] = []
    raw = Path(species_value)

    if raw.is_absolute():
        candidates.append(raw)
        if raw.suffix != ".json":
            candidates.append(raw.with_suffix(".json"))
    else:
        for prefix in (config_dir, config_dir / "ForceField"):
            candidates.append(prefix / raw)
            if raw.suffix != ".json":
                candidates.append(prefix / raw.with_suffix(".json"))

    for candidate in candidates:
        if candidate.exists():
            return candidate

    tried = "\n  ".join(str(candidate) for candidate in candidates)
    raise InputError(f"Could not find adsorbate species file {species_value!r}. Tried:\n  {tried}")


def load_config(config_path: Path) -> dict:
    try:
        with config_path.open("rb") as handle:
            return tomllib.load(handle)
    except FileNotFoundError as exc:
        raise InputError(f"Config file not found: {config_path}") from exc
    except tomllib.TOMLDecodeError as exc:
        raise InputError(f"Invalid TOML config {config_path}: {exc}") from exc


def require_section(config: dict, name: str) -> dict:
    section = config.get(name)
    if not isinstance(section, dict):
        raise InputError(f"Missing required [{name}] section in config")
    return section


def require_value(mapping: dict, key: str, context: str):
    if key not in mapping:
        raise InputError(f"Missing required config value: {context}.{key}")
    return mapping[key]


def parse_cif(cif_path: Path) -> CifData:
    lines = cif_path.read_text().splitlines()
    cell_values: dict[str, float] = {}

    for line in lines:
        parts = shlex.split(line, comments=False, posix=True)
        if len(parts) >= 2 and parts[0] in CELL_KEYS:
            cell_values[CELL_KEYS[parts[0]]] = parse_number(parts[1])

    missing = [name for name in CELL_KEYS.values() if name not in cell_values]
    if missing:
        raise InputError(f"CIF {cif_path} is missing cell parameters: {', '.join(missing)}")

    atoms = parse_atom_site_loop(lines, cif_path)
    if not atoms:
        raise InputError(f"CIF {cif_path} does not contain any atom-site rows")

    return CifData(atoms=atoms, **cell_values)


def parse_atom_site_loop(lines: list[str], cif_path: Path) -> list[Atom]:
    i = 0
    while i < len(lines):
        if lines[i].strip() != "loop_":
            i += 1
            continue

        i += 1
        headers: list[str] = []
        while i < len(lines):
            stripped = lines[i].strip()
            if not stripped:
                i += 1
                continue
            if stripped.startswith("_"):
                headers.append(stripped.split()[0])
                i += 1
                continue
            break

        if not any(header.startswith("_atom_site_") for header in headers):
            continue

        header_index = {header: idx for idx, header in enumerate(headers)}
        required = [
            "_atom_site_fract_x",
            "_atom_site_fract_y",
            "_atom_site_fract_z",
            "_atom_site_charge",
        ]
        missing = [header for header in required if header not in header_index]
        if missing:
            raise InputError(f"CIF atom-site loop in {cif_path} is missing: {', '.join(missing)}")

        name_header = None
        for candidate in ("_atom_site_type_symbol", "_atom_site_label"):
            if candidate in header_index:
                name_header = candidate
                break
        if name_header is None:
            raise InputError(f"CIF atom-site loop in {cif_path} is missing atom type/label")

        atoms: list[Atom] = []
        while i < len(lines):
            stripped = lines[i].strip()
            if not stripped or stripped == "loop_" or stripped.startswith("_") or stripped.startswith("data_"):
                break
            if stripped.startswith("#"):
                i += 1
                continue

            values = shlex.split(stripped, comments=False, posix=True)
            if len(values) < len(headers):
                raise InputError(f"Malformed atom-site row in {cif_path}: {stripped}")

            try:
                raw_name = values[header_index[name_header]]
                atoms.append(
                    Atom(
                        name=normalize_atom_name(raw_name),
                        frac_x=parse_number(values[header_index["_atom_site_fract_x"]]),
                        frac_y=parse_number(values[header_index["_atom_site_fract_y"]]),
                        frac_z=parse_number(values[header_index["_atom_site_fract_z"]]),
                        charge=parse_number(values[header_index["_atom_site_charge"]]),
                    )
                )
            except (ValueError, IndexError) as exc:
                raise InputError(f"Could not parse atom-site row in {cif_path}: {stripped}") from exc
            i += 1

        return atoms

    raise InputError(f"No atom-site loop found in CIF {cif_path}")


def normalize_atom_name(value: str) -> str:
    value = value.strip().strip("'\"")
    match = re.match(r"([A-Z][a-z]?)", value)
    if not match:
        raise ValueError(f"cannot determine atom name from {value!r}")
    return match.group(1)


def load_force_field(force_field_path: Path) -> tuple[dict[str, float], dict[str, float], dict[str, tuple[float, float]]]:
    try:
        data = json.loads(force_field_path.read_text())
    except FileNotFoundError as exc:
        raise InputError(f"Force-field file not found: {force_field_path}") from exc
    except json.JSONDecodeError as exc:
        raise InputError(f"Invalid JSON force-field file {force_field_path}: {exc}") from exc

    masses: dict[str, float] = {}
    charges: dict[str, float] = {}
    for item in data.get("PseudoAtoms", []):
        name = item.get("name")
        if not name:
            continue
        if "mass" in item:
            masses[name] = float(item["mass"])
        if "charge" in item:
            charges[name] = float(item["charge"])

    lj_params: dict[str, tuple[float, float]] = {}
    for item in data.get("SelfInteractions", []):
        name = item.get("name")
        params = item.get("parameters", item.get("params"))
        if name and params and len(params) >= 2:
            epsilon = float(params[0])
            sigma = float(params[1])
            lj_params[name] = (epsilon, sigma)

    return masses, charges, lj_params


def load_adsorbate_species(species_path: Path) -> list[AdsorbateSite]:
    try:
        data = json.loads(species_path.read_text())
    except FileNotFoundError as exc:
        raise InputError(f"Adsorbate species file not found: {species_path}") from exc
    except json.JSONDecodeError as exc:
        raise InputError(f"Invalid adsorbate species JSON {species_path}: {exc}") from exc

    pseudo_atoms = data.get("pseudoAtoms", data.get("PseudoAtoms"))
    if not isinstance(pseudo_atoms, list) or not pseudo_atoms:
        raise InputError(f"Adsorbate species file {species_path} does not contain pseudoAtoms")

    sites: list[AdsorbateSite] = []
    for item in pseudo_atoms:
        if not isinstance(item, list) or len(item) != 2:
            raise InputError(f"Malformed pseudoAtoms entry in {species_path}: {item!r}")
        name, coords = item
        if not isinstance(name, str) or not isinstance(coords, list) or len(coords) != 3:
            raise InputError(f"Malformed pseudoAtoms entry in {species_path}: {item!r}")
        try:
            sites.append(
                AdsorbateSite(
                    name=name,
                    x=float(coords[0]),
                    y=float(coords[1]),
                    z=float(coords[2]),
                )
            )
        except (TypeError, ValueError) as exc:
            raise InputError(f"Invalid pseudoAtom coordinates in {species_path}: {item!r}") from exc

    return sites


def validate_atom_types(atoms: list[Atom], masses: dict[str, float], lj_params: dict[str, tuple[float, float]]) -> None:
    atom_types = sorted({atom.name for atom in atoms})
    missing_mass = [name for name in atom_types if name not in masses]
    missing_lj = [name for name in atom_types if name not in lj_params]
    if missing_mass:
        raise InputError(f"Force field is missing PseudoAtoms mass for: {', '.join(missing_mass)}")
    if missing_lj:
        raise InputError(f"Force field is missing SelfInteractions LJ params for: {', '.join(missing_lj)}")


def validate_adsorbate_sites(
    sites: list[AdsorbateSite],
    masses: dict[str, float],
    charges: dict[str, float],
    lj_params: dict[str, tuple[float, float]],
) -> None:
    site_types = sorted({site.name for site in sites})
    missing_mass = [name for name in site_types if name not in masses]
    missing_charge = [name for name in site_types if name not in charges]
    missing_lj = [name for name in site_types if name not in lj_params]
    if missing_mass:
        raise InputError(f"Force field is missing PseudoAtoms mass for adsorbate site(s): {', '.join(missing_mass)}")
    if missing_charge:
        raise InputError(f"Force field is missing PseudoAtoms charge for adsorbate site(s): {', '.join(missing_charge)}")
    if missing_lj:
        raise InputError(f"Force field is missing SelfInteractions LJ params for adsorbate site(s): {', '.join(missing_lj)}")


def build_adsorbate_block(
    sites: list[AdsorbateSite],
    masses: dict[str, float],
    charges: dict[str, float],
    lj_params: dict[str, tuple[float, float]],
) -> tuple[list[str], float]:
    lines = [
        "------------------Adsorbate------------------",
        "Number of sites",
        str(len(sites)),
        "x(A)\ty(A)\tz(A)\tEpsilon(K)\tSigma(A) Charge(e)\tMass(g/mol)",
    ]
    total_mass = 0.0
    for site in sites:
        epsilon, sigma = lj_params[site.name]
        mass = masses[site.name]
        charge = charges[site.name]
        total_mass += mass
        lines.append(
            " ".join(
                [
                    fmt(site.x),
                    fmt(site.y),
                    fmt(site.z),
                    fmt(epsilon),
                    fmt(sigma),
                    fmt(charge),
                    fmt(mass),
                ]
            )
        )
    return lines, total_mass


def build_input(
    cif: CifData,
    atoms: list[Atom],
    masses: dict[str, float],
    lj_params: dict[str, tuple[float, float]],
    adsorbate_block: list[str],
    adsorbate_mass: float,
    config: dict,
) -> str:
    run = require_section(config, "run")
    string = require_section(config, "string")
    legacy = require_section(config, "legacy")

    lines: list[str] = [
        "Nmaxa Nmaxb Nmaxc:",
        f"{int(require_value(legacy, 'nmaxa', 'legacy'))} {int(require_value(legacy, 'nmaxb', 'legacy'))} {int(require_value(legacy, 'nmaxc', 'legacy'))}",
        "La Lb Lc dL",
        f"{fmt(cif.cell_a)} {fmt(cif.cell_b)} {fmt(cif.cell_c)} {fmt(float(require_value(legacy, 'dL', 'legacy')))}",
        "Alpha Beta Gamma",
        f"{fmt(cif.alpha)}\t{fmt(cif.beta)}\t{fmt(cif.gamma)}",
        "cutoff(A) FH_signal Mass(g/mol) Tempearture(K)  Running_steps",
        (
            f"{fmt(float(require_value(run, 'cutoff', 'run')))}\t"
            f"{int(require_value(run, 'fh_signal', 'run'))}\t"
            f"{fmt(adsorbate_mass)}\t"
            f"{fmt(float(require_value(run, 'temperature', 'run')))}\t"
            f"{int(require_value(run, 'running_steps', 'run'))}"
        ),
        "---------String Calculation Settings---------",
        "Direction",
        str(int(require_value(string, "direction", "string"))),
        "#_of_points delta_frac delta_angle_degree",
        (
            f"{int(require_value(string, 'points', 'string'))} "
            f"{fmt(float(require_value(string, 'delta_frac', 'string')))} "
            f"{fmt(float(require_value(string, 'delta_angle_degree', 'string')))}"
        ),
        "convergence_setting",
        str(require_value(string, "convergence_setting", "string")),
    ]

    lines.extend(adsorbate_block)
    lines.extend(
        [
            "------------------Adsorbent------------------",
            "Number of atoms",
            str(len(atoms)),
            "ID diameter(A) Epsilon(K) Charge(e) mass(g/mol) frac_x frac_y frac_z atom_name",
        ]
    )

    for index, atom in enumerate(atoms, start=1):
        epsilon, sigma = lj_params[atom.name]
        lines.append(
            "\t".join(
                [
                    str(index),
                    fmt(sigma),
                    fmt(epsilon),
                    fmt(atom.charge),
                    fmt(masses[atom.name]),
                    fmt(atom.frac_x),
                    fmt(atom.frac_y),
                    fmt(atom.frac_z),
                    atom.name,
                ]
            )
        )

    return "\n".join(lines) + "\n"


def fmt(value: float) -> str:
    return f"{value:.10g}"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Generate a GaSSM sim.input file from a CIF and force-field config.")
    parser.add_argument("config", help="Path to TOML config file, for example sim_config.toml")
    args = parser.parse_args(argv)

    try:
        config_path = Path(args.config).resolve()
        config = load_config(config_path)
        config_dir = config_path.parent

        structure_name = str(require_value(config, "structure", "config"))
        output_path = resolve_path(str(require_value(config, "output", "config")), config_dir)
        force_field_path = resolve_path(str(config.get("force_field", "ForceField/force_field.json")), config_dir)
        adsorbate = require_section(config, "adsorbate")
        adsorbate_species = str(require_value(adsorbate, "species", "adsorbate"))
        cif_path = resolve_structure(structure_name, config_dir)
        species_path = resolve_adsorbate_species(adsorbate_species, config_dir)

        cif = parse_cif(cif_path)
        masses, charges, lj_params = load_force_field(force_field_path)
        validate_atom_types(cif.atoms, masses, lj_params)
        adsorbate_sites = load_adsorbate_species(species_path)
        validate_adsorbate_sites(adsorbate_sites, masses, charges, lj_params)
        adsorbate_block, adsorbate_mass = build_adsorbate_block(adsorbate_sites, masses, charges, lj_params)

        output = build_input(cif, cif.atoms, masses, lj_params, adsorbate_block, adsorbate_mass, config)
        output_path.write_text(output)

        atom_types = ", ".join(sorted({atom.name for atom in cif.atoms}))
        direction = require_section(config, "string").get("direction")
        print(f"Wrote {output_path}")
        print(f"Structure: {cif_path.name}")
        print(f"Adsorbate: {species_path.stem}")
        print(f"Framework atoms: {len(cif.atoms)}")
        print(f"Atom types: {atom_types}")
        print(f"Direction: {direction}")
        return 0
    except InputError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
