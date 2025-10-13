#!/usr/bin/env python3

import numpy as np
from os.path import join
import oer
import yaml
import re

MOL_SPECIES = {"NO3RR": ["H2", "H2O", "H3N", "HNO3"]}

def path_map(yaml_file="symlinks.yaml"):
    with open(yaml_file, "r") as f:
        data = yaml.safe_load(f)
    return data


def mongo_composition_match(field, composition):
    """
    Generate a list of MongoDB query clauses for exact, order-neutral matching
    of a composition dictionary stored at `field` (e.g. "output.composition").

    Args:
        field (str): The field path to the composition dict in MongoDB documents.
        composition (dict): The target composition, e.g. {"C": 1, "O": 2}

    Returns:
        list: List of query components (to be combined with $and)
    """
    queries = [{f"{field}.{k}": v for k, v in composition.items()}]
    queries.append(
        {
            "$expr": {
                "$eq": [{"$size": {"$objectToArray": f"${field}"}}, len(composition)]
            }
        }
    )
    # Flatten the dict for element queries
    element_queries = [{f"{field}.{k}": v} for k, v in composition.items()]
    return element_queries + [queries[-1]]


def get_adsorbate_energy(
    store,
    functional,
    series,
    aexx=None,
    thermo_corr=False
):
    gga, lvdw, lhfcalc = oer.get_param(functional)
    query = [
        {"name": "adsorbate relax"},
        {"output.input.parameters.GGA": {"$regex": gga, "$options": "i"}},
        {"output.input.parameters.LUSE_VDW": lvdw},
        {"output.input.parameters.LHFCALC": lhfcalc},
    ]
    if lhfcalc and aexx:
        query.append({"output.input.parameters.AEXX": aexx})
    if isinstance(series, str):
        query += [{"output.dir_name": regex_dir_name(series)}]
    elif isinstance(series, list):
        for _s in series:
            query += [{"output.dir_name": regex_dir_name(_s)}]
    doc = store.query_one({"$and": query})
    return doc["output"]["output"]["energy"] + doc["thermo_corr"] if thermo_corr else doc["output"]["output"]["energy"]


def report_noxrr(
    store,
    substrate_string,
    functional,
    reaction="NO3RR",
    react_coords=None,
    n_desorbed=None,
    n_protons=None,
    series=None,
    series_mol=None,
    aexx=None,
    thermo_corr=False,
    write_yaml=True,
    yaml_prefix="noxrr",
):
    gga, lvdw, lhfcalc = oer.get_param(functional)
    energy = {}
    struct = {}
    forces = {}
    for i in MOL_SPECIES[reaction]:
        energy[i], struct[i], forces[i], _ = oer.get_energy_and_structure(
            store, i, functional, series_mol, aexx, thermo_corr
        )
    energy["H"] = 0.5 * energy["H2"]

    # substrate
    energy["substrate"], struct["substrate"], forces["substrate"], comp_substrate, _ = (
        oer.find_substrate(store, substrate_string, functional, series, aexx, thermo_corr)
    )

    # go through along the reaction coords
    delta_g = {}
    free_energy = {}
    e0 = energy["substrate"] + energy["HNO3"]

    for i, rc in enumerate(react_coords):
        rc_path = ([series] if isinstance(series, str) else list(series)) + [rc]
        energy_rc = get_adsorbate_energy(store, functional, rc_path, thermo_corr=thermo_corr)
        #print(rc, energy_rc)
        free_energy["*" + rc] = energy_rc
        delta_g["*" + rc] = (energy_rc
            + energy["H2O"] * n_desorbed[i]
            - energy["H"] * n_protons[i] - e0
        )

    print(free_energy, delta_g)

    if write_yaml:
        _to_yaml(
            yaml_prefix=yaml_prefix,
            react_coords=react_coords,
            n_desorbed=n_desorbed,
            n_protons=n_protons,
            free_energy=free_energy,
            delta_g=delta_g,
            reaction=reaction,
        )


def get_energies(doc):
    return (
        doc["output"]["output"]["output"]["final_energy"],
        doc["output"]["output"]["output"]["free_energy"],
    )


def regex_dir_name(name):
    pattern = rf"(^|/){re.escape(name)}(/|$)"
    return re.compile(pattern)


def report_noxrr_mace(
    store,
    substrate_string,
    reaction="NO3RR",
    react_coords=None,
    n_desorbed=None,
    n_protons=None,
    series=None,
    write_yaml=True,
    yaml_prefix="noxrr",
):
    from pymatgen.core import Composition

    energy = {}
    free_energy = {}

    comp_substrate = Composition(substrate_string).as_dict()
    query = [{"name": "substrate relax"}, {"output.input.calculator": "macecalculator"}]
    query += mongo_composition_match("output.composition", comp_substrate)
    if isinstance(series, str):
        query += [{"output.dir_name": regex_dir_name(series)}]
    elif isinstance(series, list):
        for _s in series:
            query += [{"output.dir_name": regex_dir_name(_s)}]
    docs = [i for i in store.query({"$and": query})]
    if len(docs) > 1:
        raise ValueError("More than one substrate entry found!")
    energy["substrate"], free_energy["substrate"] = get_energies(docs[0])

    # molecules
    for mol in MOL_SPECIES[reaction]:
        comp_mol = Composition(mol).as_dict()
        query = [
            {"name": "molecule relax"},
            {"output.input.calculator": "macecalculator"},
        ]
        query += mongo_composition_match("output.composition", comp_mol)
        docs = [i for i in store.query({"$and": query})]
        try:
            energy[mol], free_energy[mol] = get_energies(docs[0])
        except Exception as e:
            print(f"{mol} not found")

    delta_g = {}
    # the first step is always the adsorption of CO2
    for i, rc in enumerate(react_coords):
        query = [
            {"name": "adsorbate relax"},
            {"output.input.calculator": "macecalculator"},
        ]
        query += [{"output.dir_name": regex_dir_name(rc)}]
        if isinstance(series, str):
            query += [{"output.dir_name": regex_dir_name(series)}]
        elif isinstance(series, list):
            for _s in series:
                query += [{"output.dir_name": regex_dir_name(_s)}]
        docs = [i for i in store.query({"$and": query})]
        if len(docs) > 1:
            raise ValueError("More than one adsorbate entry found!")
        energy["*" + rc], free_energy["*" + rc] = get_energies(docs[0])

    for i, rc in enumerate(react_coords):
        delta_g["*" + rc] = (
            free_energy["*" + rc]
            - free_energy["HNO3"]
            - free_energy["substrate"]
            + free_energy["H2O"] * n_desorbed[i]
            - 0.5 * free_energy["H2"] * n_protons[i]
        )
    last = react_coords[-1].split('_')[0]
    if last.upper() in MOL_SPECIES[reaction]:
        delta_g[last + "_g"] = (
            free_energy[last.upper()]
            + free_energy["H2O"] * n_desorbed[-1]
            - 0.5 * free_energy["H2"] * n_protons[-1]
            - free_energy["HNO3"]
        )

    if write_yaml:
        _to_yaml(
            yaml_prefix=yaml_prefix,
            react_coords=react_coords,
            n_desorbed=n_desorbed,
            n_protons=n_protons,
            energy=energy,
            free_energy=free_energy,
            delta_g=delta_g,
            reaction=reaction,
        )

    return energy, free_energy, delta_g


def _to_yaml(
    yaml_prefix, react_coords, n_desorbed, n_protons, free_energy, delta_g, reaction
):
    data = {
        rc: {
            "n_desorbed": nd,
            "n_protons": np,
            "free_energy": free_energy["*" + rc],
            "delta_g": delta_g["*" + rc],
        }
        for rc, nd, np in zip(react_coords, n_desorbed, n_protons)
    }
    last = react_coords[-1].split('_')[0]
    if last.upper() in MOL_SPECIES[reaction]:
        data[last+'_g'] = {
                "n_desorbed": n_desorbed[-1],
                "n_protons": n_protons[-1],
                "delta_g": delta_g[last+'_g'],
            }
        
    with open(f"{yaml_prefix}.yaml", "w") as f:
        yaml.dump(data, f, sort_keys=False)


if __name__ == "__main__":
    store = oer.connect_db()
    substrate = "In36 Ga12 Cu48 S96" 
    react_coords = ["hno3_3", "hno3-h_3", "no2_3", "no2-h_3", "no2-h-h_3", "no_3", "no-h_3", "no-h-h_3", "nh_oh-h_3", "n-h_3", "n-h-h_3", "n-h-h-h_3"]
    n_desorbed = [0, 0, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3]
    n_protons =  [0, 1, 1, 2, 3, 3, 4, 5, 6, 6, 7, 8]

    energy, free_energy, delta_g = report_noxrr_mace(
        store,
        substrate,
        series="cigs_112_cat",
        react_coords=react_coords,
        n_desorbed=n_desorbed,
        n_protons=n_protons,
    )
