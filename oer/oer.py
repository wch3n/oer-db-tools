#!/usr/bin/env python3

import numpy as np
from os.path import join
from jobflow import SETTINGS
from pymatgen.core import Structure, Composition, Molecule
import re

CORR = {}
CORR["H2"] = 0.268 - 0.554 + 0.039 + 0.026 + 0.026
CORR["O2"] = 0.097 - 0.730 + 0.039 + 0.026 + 0.026
CORR["H2O"] = 0.576 - 0.720 + 0.039 + 0.026 + 0.018 + 0.026
CORR["H2O"] -= 0.20 # gas to liquid
CORR["OH*"] = 0.350 - 0.086 + 0.051 
CORR["O*"] = 0.072 - 0.078 + 0.038 
CORR["OOH*"] = 0.467 - 0.158 + 0.081 
CORR["H*"] = 0.175 - 0.015 + 0.011
ADSORBATE_SPECIES = {"OER": ["O", "OH", "OOH"], "OER_bi": ["O", "OH", "H"], "HER": ["H"]}
DIS_TOL_MAX = 1.0


def connect_db():
    store = SETTINGS.JOB_STORE
    store.connect()
    return store


def get_energy_and_structure(store, formula, functional="PBE", series=None):
    gga, lvdw, lhfcalc = get_param(functional)
    query =[
                {"output.formula_pretty": formula},
                {"output.input.parameters.LUSE_VDW": lvdw},
                {"output.input.parameters.GGA": gga},
                {"output.input.parameters.LHFCALC": lhfcalc},
            ]
    if isinstance(series, str):
        query.append({"output.dir_name": {"$regex": f'/{series}/'}})
    doc = store.query_one(
        {
            "$and": query
        }
    )
    print(formula, doc["output"]["output"]["energy"], doc["output"]["dir_name"])
    return (
        doc["output"]["output"]["energy"],
        Structure.from_dict(doc["output"]["structure"]),
        doc["output"]["output"]["forces"],
        doc,
    )


def calc_deltaG(energy, reaction="OER", corr=True):
    if corr:
        energy = apply_corr(energy)
    match reaction:
        case "OER":
            g1 = (
                energy["OH*"] - energy["substrate"] - energy["H2O"] + 0.5 * energy["H2"]
            )
            g2 = energy["O*"] - energy["OH*"] + 0.5 * energy["H2"]
            g3 = energy["OOH*"] - energy["O*"] - energy["H2O"] + 0.5 * energy["H2"]
            g4 = (
                energy["substrate"] + energy["O2"] - energy["OOH*"] + 0.5 * energy["H2"]
            )
            return {"G1": g1, "G2": g2, "G3": g3, "G4": g4}
        case "OER_bi":
            g1 = (
                energy["OH*"] - energy["substrate"] - energy["H2O"] + 0.5 * energy["H2"]
            )
            g2 = energy["O*"] - energy["OH*"] + 0.5 * energy["H2"]
            g3 = energy["H*"] + energy["O2"] - energy["O*"] - energy["H2O"] + 0.5 * energy["H2"]
            g4 = (
                energy["substrate"] - energy["H*"] + 0.5 * energy["H2"]
            )
            return {"G1": g1, "G2": g2, "G3": g3, "G4": g4}
        case "HER":
            g1 = energy["H*"] - energy["substrate"] - 0.5 * energy["H2"]
            return {"G1": g1}


def get_param(functional):
    if functional == "PBE":
        gga = "PE"
        lvdw = False
        lhfcalc = False
    elif functional == "rVV10":
        gga = "ML"
        lvdw = True
        lhfcalc = False
    elif functional == "PBE0":
        gga = "PE"
        lvdw = False
        lhfcalc = True
    return gga, lvdw, lhfcalc


def report(
    store,
    substrate_string,
    functional,
    adsorbate_index=None,
    anchor_index=None,
    write_poscar=False,
    exhaustive=True,
    reaction="OER",
    series=None,
    series_mol = None,
):
    gga, lvdw, lhfcalc = get_param(functional)
    energy = {}
    struct = {}
    forces = {}
    for i in ["O2", "H2", "H2O"]:
        energy[i], struct[i], forces[i], _ = get_energy_and_structure(
            store, i, functional, series_mol
        )

    # substrate
    energy["substrate"], struct["substrate"], forces["substrate"], comp_substrate, _ = (
        find_substrate(store, substrate_string, functional, series)
    )

    # adsorbate
    query = [
        {"name": "adsorbate relax"},
        {"output.input.parameters.GGA": gga},
        {"output.input.parameters.LUSE_VDW": lvdw},
        {"output.input.parameters.LHFCALC": lhfcalc},
    ]
    if isinstance(series, str):
        #query.append({"output.dir_name": {"$regex": re.escape(series)}})
        query.append({"output.dir_name": {"$regex": f'/{series}/'}})
    for ads in ADSORBATE_SPECIES[reaction]:
        tot_energy = 0.0
        comp = comp_substrate + Composition(ads)
        for adsorbate in store.query({"$and": query}):
            comp_adsorbate = Composition(adsorbate["output"]["composition"])
            if comp_adsorbate.reduced_composition == comp.reduced_composition:
                s = Structure.from_dict(adsorbate["output"]["structure"])
                ads_indices, is_wrong_adsorbate = find_adsorbate(
                    s, struct["substrate"], ads
                )
                if is_wrong_adsorbate:
                    # print('Adsorbate different from the specified species...')
                    continue
                if anchor_index is not None:
                    binding_site, adsorbate_site, binding_dist = find_anchor(
                        s, ads_indices
                    )
                    if (
                        isinstance(anchor_index, int) and binding_site == anchor_index
                    ) or (
                        isinstance(anchor_index, list) and binding_site in anchor_index
                    ):
                        # print(ads_indices, binding_dist)
                        print(
                            ads,
                            binding_site,
                            "-".join(
                                [
                                    s[binding_site].species_string,
                                    s[adsorbate_site].species_string,
                                ]
                            ),
                            adsorbate["output"]["output"]["energy"],
                            f"{adsorbate['output']['dir_name']}",
                        )
                    else:
                        # print('Wrong adsorption site...')
                        continue
                else:
                    raise NotImplementedError

                if adsorbate["output"]["output"]["energy"] < tot_energy:
                    (
                        energy["".join([ads, "*"])],
                        struct["".join([ads, "*"])],
                        forces["".join([ads, "*"])],
                    ) = (
                        adsorbate["output"]["output"]["energy"],
                        Structure.from_dict(adsorbate["output"]["structure"]),
                        adsorbate["output"]["output"]["forces"],
                    )
                    tot_energy = adsorbate["output"]["output"]["energy"]
                if not exhaustive:
                    break

    deltaG = calc_deltaG(energy, reaction, corr=True)
    if write_poscar:
        save_poscar(struct)
    return deltaG, energy, struct, forces


def apply_corr(energy):
    for i in energy.keys():
        if i in CORR.keys():
            energy[i] += CORR[i]
    return energy


def extract_elements(formula):
    pattern = r"[A-Z][a-z]?"
    elements = re.findall(pattern, formula)
    return elements


def find_adsorbate(ads, substrate_structure, adsorbate_string):
    elements = extract_elements(adsorbate_string)
    # increase dis_tol up to DIS_TOL_MAX if adsorbate not found
    while dis_tol < DIS_TOL_MAX:
        print(dis_tol)
        subs_matched = []
        for i in range(len(ads)):
            for j in substrate_structure:
                if ads[i].distance(j) < dis_tol:
                    subs_matched.append(i)
                    break
        ads_indices = set(range(len(ads))).symmetric_difference(subs_matched)
        ads_indices = [i for i in ads_indices if ads[i].species_string in elements]
        mol = Molecule.from_sites([ads[i] for i in ads_indices])
        comp_ads = Composition(adsorbate_string)
        if mol.composition == comp_ads:
            break
        else:
            dis_tol += 0.2
    bonds = ["".join(sorted(elements[i : i + 2])) for i in range(len(elements) - 1)]
    bonds_found = [
        "".join(sorted([i.site1.species_string, i.site2.species_string]))
        for i in mol.get_covalent_bonds()
    ]
    is_wrong_adsorbate = True if set(bonds_found) - set(bonds) else False
    return ads_indices, is_wrong_adsorbate


def write_poscar(struct, adsorbate_only=True):
    for name in struct:
        if adsorbate_only and not "*" in name:
            continue
        struct[name].to_file(f'POSCAR.{name.replace("*","-")}')


def find_substrate(store, substrate_string, functional, series=None):
    gga, lvdw, lhfcalc = get_param(functional)
    comp_target = Composition(substrate_string)
    query = [
        {"name": "substrate relax"},
        {"output.input.parameters.GGA": gga},
        {"output.input.parameters.LUSE_VDW": lvdw},
        {"output.input.parameters.LHFCALC": lhfcalc},
    ]
    if isinstance(series, str):
        #query.append({"output.dir_name": {"$regex": re.escape(series)}})
        query.append({"output.dir_name": {"$regex": f'/{series}/'}})
    for substrate in store.query({"$and": query}):
        comp_substrate = Composition(substrate["output"]["composition"])
        if comp_substrate.reduced_composition == comp_target.reduced_composition:
            # print(f"{'SLAB':<6}", f"{-1:<4}", f"{'XX-XX':>2}", ' '.join(f'{x:8.2f}' for x in [0,0]),
            #    f'{substrate["output"]["output"]["energy"]:8.2f}', f"{substrate['output']['dir_name']}")
            energy, struct, forces, doc = get_energy_and_structure(
                store, substrate["output"]["formula_pretty"], functional, series
            )
            break
    return energy, struct, forces, comp_substrate, doc


def find_anchor(structure, ads_indices):
    binding_dist = 10.0
    binding_site = None
    adsorbate_site = None
    for idx in ads_indices:
        _d = structure.distance_matrix[idx]
        mask = np.zeros(_d.shape, dtype=bool)
        mask[ads_indices] = True
        _d_masked = np.ma.masked_array(_d, mask=mask)
        _site = np.argsort(_d_masked)[0]
        if _d[_site] < binding_dist:
            binding_dist = _d[_site]
            binding_site = _site
            adsorbate_site = idx
    return binding_site, adsorbate_site, binding_dist


def find_all(store, substrate_string, functional, reaction="OER", series=None, series_mol=None):
    gga, lvdw, lhfcalc = get_param(functional)
    energy = {}
    struct = {}
    forces = {}
    docs = {}
    for i in ["O2", "H2", "H2O"]:
        energy[i], struct[i], forces[i], docs[i] = get_energy_and_structure(
            store, i, functional, series_mol
        )
    comp_target = Composition(substrate_string)
    (
        energy["substrate"],
        struct["substrate"],
        forces["substrate"],
        comp_substrate,
        docs["SLAB"],
    ) = find_substrate(store, substrate_string, functional, series)
    for ads in ADSORBATE_SPECIES[reaction]:
        comp = comp_substrate + Composition(ads)
        query = [
            {"name": "adsorbate relax"},
            {"output.input.parameters.GGA": gga},
            {"output.input.parameters.LUSE_VDW": lvdw},
            {"output.input.parameters.LHFCALC": lhfcalc},
        ]
        if isinstance(series, str):
            #query.append({"output.dir_name": {"$regex": re.escape(series)}})
            query.append({"output.dir_name": {"$regex": f'/{series}/'}})
        for adsorbate in store.query({"$and": query}):
            comp_adsorbate = Composition(adsorbate["output"]["composition"])
            if comp_adsorbate.reduced_composition == comp.reduced_composition:
                s = Structure.from_dict(adsorbate["output"]["structure"])
                ads_indices, is_wrong_adsorbate = find_adsorbate(
                    s, struct["substrate"], ads
                )
                if is_wrong_adsorbate:
                    # print('Adsorbate different from the specified species...')
                    continue
                binding_site, adsorbate_site, _ = find_anchor(s, ads_indices)
                energy["".join([ads, "*"])], struct["".join([ads, "*"])] = adsorbate[
                    "output"
                ]["output"]["energy"], Structure.from_dict(
                    adsorbate["output"]["structure"]
                )
                pair = f"{s[binding_site].species_string:>2}-{s[adsorbate_site].species_string:<2}"
                ads_name = f"{ads+'*':<6}"
                total_energy = f'{adsorbate["output"]["output"]["energy"]:8.2f}'
                print(
                    ads_name,
                    f"{binding_site:<4}",
                    pair,
                    " ".join(f"{x:8.2f}" for x in s[adsorbate_site].coords[:2]),
                    total_energy,
                    f"{adsorbate['output']['dir_name']}",
                )
                docs[ads + "*"] = adsorbate
    return docs


def save_poscar(struct, adsorbate_only=True):
    for name in struct:
        if adsorbate_only and not "*" in name:
            continue
        struct[name].to_file(f'POSCAR.{name.replace("*","-")}')
