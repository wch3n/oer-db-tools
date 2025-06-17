#!/usr/bin/env python3

import numpy as np
from os.path import join
import oer
import yaml
import re

MOL_SPECIES = {"CO2RR": ["CO2", "H2", "H2O", "CO", "CH3OH", "CH4", "HCOOH"]}


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


def find_adsorbate(
    store,
    functional,
    series,
    aexx=None,
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
    return doc["output"]["output"]["energy"]


def report_co2rr(
    store,
    substrate_string,
    functional,
    reaction="CO2RR",
    react_coords=None,
    n_desorbed=None,
    n_protons=None,
    series=None,
    series_mol=None,
    aexx=None,
):
    gga, lvdw, lhfcalc = oer.get_param(functional)
    energy = {}
    struct = {}
    forces = {}
    for i in MOL_SPECIES[reaction]:
        energy[i], struct[i], forces[i], _ = oer.get_energy_and_structure(
            store, i, functional, series_mol, aexx
        )
    energy["H"] = 0.5 * energy["H2"]

    # substrate
    energy["substrate"], struct["substrate"], forces["substrate"], comp_substrate, _ = (
        oer.find_substrate(store, substrate_string, functional, series, aexx)
    )

    # go through along the reaction coords
    react_coords_path = path_map()
    delta_g = []
    delta_g_raw = []
    delta_g.append(energy["substrate"] + energy["CO2"])
    delta_g_raw.append(energy["substrate"] + energy["CO2"])

    for i, rc in enumerate(react_coords[1:]):
        print(rc)
        rc_path = react_coords_path[rc]
        delta_g_raw.append(find_adsorbate(store, functional, series + [rc_path]))
        delta_g.append(
            find_adsorbate(store, functional, series + [rc_path])
            + energy["H2O"] * n_desorbed[i + 1]
            - energy["H"] * n_protons[i + 1]
        )

    print([i for i in zip(react_coords, np.array(delta_g_raw))])
    print([i for i in zip(react_coords, np.array(delta_g) - delta_g[0])])


def get_energies(doc):
    return (
        doc["output"]["output"]["output"]["final_energy"],
        doc["output"]["output"]["output"]["free_energy"],
    )


def regex_dir_name(name):
    pattern = rf"(^|/){re.escape(name)}(/|$)"
    return re.compile(pattern)


def report_co2rr_mace(
    store,
    substrate_string,
    reaction="CO2RR",
    react_coords=None,
    n_desorbed=None,
    n_protons=None,
    series=None,
    write_yaml=True,
    yaml_prefix="co2rr",
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
            - free_energy["CO2"]
            - free_energy["substrate"]
            + free_energy["H2O"] * n_desorbed[i]
            - 0.5 * free_energy["H2"] * n_protons[i]
        )
    if react_coords[-1].upper() in MOL_SPECIES[reaction]:
        delta_g[react_coords[-1] + "_g"] = (
            free_energy[react_coords[-1].upper()]
            + free_energy["H2O"] * n_desorbed[-1]
            - 0.5 * free_energy["H2"] * n_protons[-1]
            - free_energy["CO2"]
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
    yaml_prefix, react_coords, n_desorbed, n_protons, energy, free_energy, delta_g, reaction
):
    data = {
        rc: {
            "n_desorbed": nd,
            "n_protons": np,
            "energy": energy["*" + rc],
            "free_energy": free_energy["*" + rc],
            "delta_g": delta_g["*" + rc],
        }
        for rc, nd, np in zip(react_coords, n_desorbed, n_protons)
    }
    if react_coords[-1].upper() in MOL_SPECIES[reaction]:
        data[react_coords[-1]+'_g'] = {
                "n_desorbed": n_desorbed[-1],
                "n_protons": n_protons[-1],
                "delta_g": delta_g[react_coords[-1]+'_g'],
            }
        
    with open(f"{yaml_prefix}.yaml", "w") as f:
        yaml.dump(data, f, sort_keys=False)


if __name__ == "__main__":
    store = oer.connect_db()
    substrate = "Ti90 C60 O48 H48 F12 Cu15"  # 0.33 ML
    react_coords = ["co2", "co2-h", "co2-h-h", "co", "cho", "ch2o", "ch3o", "ch3oh"]
    n_desorbed = [0, 0, 0, 1, 1, 1, 1, 1]
    n_protons = [0, 1, 2, 2, 3, 4, 5, 6]

    energy, free_energy, delta_g = report_co2rr_mace(
        store,
        substrate,
        series="ML_0_33",
        react_coords=react_coords,
        n_desorbed=n_desorbed,
        n_protons=n_protons,
    )
