#!/usr/bin/env python3

import numpy as np
from os.path import join
from jobflow import SETTINGS
from pymatgen.core import Structure, Composition, Molecule
import re

CORR = {}
CORR["H2"] = 0.268 - 0.403 + 0.039 + 0.026 + 0.026
CORR["O2"] = 0.097 - 0.634 + 0.039 + 0.026 + 0.026
CORR["H2O"] = 0.567 - 0.670 + 0.039 + 0.039 + 0.026
CORR["OH*"] = 0.350 - 0.086 + 0.051 
CORR["O*"] = 0.072 - 0.078 + 0.038 
CORR["OOH*"] = 0.467 - 0.158 + 0.081 
CORR["H*"] = 0.175 - 0.015 + 0.011
CORR["CH3OH"] = 1.343 - 0.751 + 0.039 + 0.039 + 0.024 + 0.026
ADSORBATE_SPECIES = {"OER": ["O", "OH", "OOH"], "OER_bi": ["O", "OH", "H"], "HER": ["H"],
                     "CORR_1": ["CO", "CHO", "CH2O", "CH3O", "CH3OH"],}
MOL_SPECIES = {"OER": ["H2", "O2", "H2O"], 
               "CORR_1": ["CO", "H2", "H4CO"]}
DIS_TOL_MAX = 1.0


def connect_db():
    store = SETTINGS.JOB_STORE
    store.connect()
    return store


def get_energy_and_structure(store, formula, functional="PBE", series=None, aexx=None):
    gga, lvdw, lhfcalc = get_param(functional)
    query =[
                {"output.formula_pretty": formula},
                {"output.input.parameters.LUSE_VDW": lvdw},
                {"output.input.parameters.GGA": {"$regex": gga, "$options": "i"}},
                {"output.input.parameters.LHFCALC": lhfcalc},
            ]
    if lhfcalc and aexx:
        query.append({"output.input.parameters.AEXX": aexx})
    if isinstance(series, str):
        query.append({"output.dir_name": {"$regex": f'/{series}/'}})
    elif isinstance(series, list):
        for _s in series:
            query.append({"output.dir_name": {"$regex": f'/{_s}/'}})
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


def calc_deltaG(energy_raw, reaction="OER", corr=True, corr_liquid=0):
    energy = energy_raw.copy()
    if corr:
        energy = apply_corr(energy, corr_liquid)
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
        case "CORR_1":
            g1 = energy["CHO*"] - energy["CO*"] - 0.5 * energy["H2"]
            g2 = energy["CH2O*"] - energy["CHO*"] - 0.5 * energy["H2"]
            g3 = energy["CH3O*"] - energy["CH2O*"] - 0.5 * energy["H2"]
            g4 = energy["CH3OH*"] - energy["CH3O*"] - 0.5 * energy["H2"]
            g5 = energy["CH3OH"] + substrate - energy["CH3OH*"]
            return {"G1": g1, "G2": g2, "G3": g3, "G4": g4, "G5": g5}


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
    elif functional == "PBE0-rVV10":
        gga = "PE"
        lvdw = True
        lhfcalc = True
    elif functional == "PBE-D3":
        gga = "PE"
        lvdw = False
        lhfcalc = False
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
    corr_liquid = 0,
    aexx=None,
):
    gga, lvdw, lhfcalc = get_param(functional)
    energy = {}
    struct = {}
    forces = {}
    for i in MOL_SPECIES[reaction]:
        energy[i], struct[i], forces[i], _ = get_energy_and_structure(
            store, i, functional, series_mol, aexx
        )

    # substrate
    energy["substrate"], struct["substrate"], forces["substrate"], comp_substrate, _ = (
        find_substrate(store, substrate_string, functional, series, aexx)
    )

    # adsorbate
    query = [
        {"name": "adsorbate relax"},
        {"output.input.parameters.GGA": {"$regex": gga, "$options": "i"}},
        {"output.input.parameters.LUSE_VDW": lvdw},
        {"output.input.parameters.LHFCALC": lhfcalc},
    ]
    if lhfcalc and aexx:
        query.append({"output.input.parameters.AEXX": aexx})
    if isinstance(series, str):
        query.append({"output.dir_name": {"$regex": f'/{series}/'}})
    elif isinstance(series, list):
        for _s in series:
            query.append({"output.dir_name": {"$regex": f'/{_s}/'}})
    for ads in ADSORBATE_SPECIES[reaction]:
        tot_energy = 0.0
        comp = comp_substrate + Composition(ads)
        for adsorbate in store.query({"$and": query}):
            comp_adsorbate = Composition(adsorbate["output"]["composition"])
            if comp_adsorbate == comp:
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

    deltaG = calc_deltaG(energy, reaction, corr=True, corr_liquid=corr_liquid)
    if write_poscar:
        save_poscar(struct)
    return deltaG, energy, struct, forces


def apply_corr(energy, corr_liquid):
    CORR['H2O'] += corr_liquid
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
    dis_tol = 0.2
    while dis_tol < DIS_TOL_MAX:
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


def find_substrate(store, substrate_string, functional, series=None, aexx=None):
    gga, lvdw, lhfcalc = get_param(functional)
    comp_target = Composition(substrate_string)
    query = [
        {"name": "substrate relax"},
        {"output.input.parameters.GGA": {"$regex": gga, "$options": "i"}},
        {"output.input.parameters.LUSE_VDW": lvdw},
        {"output.input.parameters.LHFCALC": lhfcalc},
    ]
    if lhfcalc and aexx:
        query.append({"output.input.parameters.AEXX": aexx})
    if isinstance(series, str):
        query.append({"output.dir_name": {"$regex": f'/{series}/'}})
    elif isinstance(series, list):
        for _s in series:
            query.append({"output.dir_name": {"$regex": f'/{_s}/'}})

    slabs = []
    for substrate in store.query({"$and": query}):
        if Composition(substrate["output"]["composition"]) == comp_target:
            # print(f"{'SLAB':<6}", f"{-1:<4}", f"{'XX-XX':>2}", ' '.join(f'{x:8.2f}' for x in [0,0]),
            #    f'{substrate["output"]["output"]["energy"]:8.2f}', f"{substrate['output']['dir_name']}")
            slabs.append(substrate)
            comp_substrate = Composition(substrate["output"]["composition"])

    eslab = 0
    for slab in slabs:
        if slab["output"]["output"]["energy"] < eslab:
            substrate = slab.copy()
            eslab = slab["output"]["output"]["energy"]

    energy, struct, forces, doc = get_energy_and_structure(
            store, substrate["output"]["formula_pretty"], functional, series, aexx
            )

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

def find_by_name(store, functional, series, name, aexx=None):
    gga, lvdw, lhfcalc = get_param(functional)
    energy = {}
    struct = {}
    forces = {}
    docs = {}

    query = [
        {"name": name},
        {"output.input.parameters.GGA": {"$regex": gga, "$options": "i"}},
        {"output.input.parameters.LUSE_VDW": lvdw},
        {"output.input.parameters.LHFCALC": lhfcalc},
    ]
    if lhfcalc and aexx:
        query.append({"output.input.parameters.AEXX": aexx})
    if isinstance(series, str):
        query.append({"output.dir_name": {"$regex": f'/{series}/'}})
    elif isinstance(series, list):
        for _s in series:
            query.append({"output.dir_name": {"$regex": f'/{_s}/'}})

    docs = []
    for entry in store.query({"$and": query}):
        docs.append(entry)
        total_energy = f'{entry["output"]["output"]["energy"]:8.2f}'
        composition = entry["output"]["composition"]
        print(composition,
            total_energy,
            f"{entry['output']['dir_name']}"
        )

    return docs

def find_all(store, substrate_string, functional, reaction="OER", series=None, series_subs=None, series_mol=None, name="adsorbate relax", aexx=None, fast_mode=False):
    gga, lvdw, lhfcalc = get_param(functional)
    energy = {}
    struct = {}
    forces = {}
    docs = {}
    for i in MOL_SPECIES[reaction]:
        energy[i], struct[i], forces[i], docs[i] = get_energy_and_structure(
            store, i, functional, series_mol, aexx
        )
    (
        energy["substrate"],
        struct["substrate"],
        forces["substrate"],
        comp_substrate,
        docs["SLAB"],
    ) = find_substrate(store, substrate_string, functional, series_subs, aexx)
        
    query = [
        {"name": name},
        {"output.input.parameters.GGA": {"$regex": gga, "$options": "i"}},
        {"output.input.parameters.LUSE_VDW": lvdw},
        {"output.input.parameters.LHFCALC": lhfcalc},
    ]
    if lhfcalc and aexx:
        query.append({"output.input.parameters.AEXX": aexx})
    if isinstance(series, str):
        query.append({"output.dir_name": {"$regex": f'/{series}/'}})
    elif isinstance(series, list):
        for _s in series:
            query.append({"output.dir_name": {"$regex": f'/{_s}/'}})

    comp = {}

    for ads in ADSORBATE_SPECIES[reaction]:
        docs[ads+'*'] = []
        comp[ads+'*'] = comp_substrate + Composition(ads)

    docs['BARE'] = []
    for adsorbate in store.query({"$and": query}):
        comp_adsorbate = Composition(adsorbate["output"]["composition"])
        for ads in ADSORBATE_SPECIES[reaction]:
            if comp_adsorbate == comp[ads+'*']:
                s = Structure.from_dict(adsorbate["output"]["structure"])
                if not fast_mode:
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
                if not fast_mode:
                    print(
                        ads_name,
                        f"{binding_site:<4}",
                        pair,
                        " ".join(f"{x:8.2f}" for x in s[adsorbate_site].coords[:2]),
                        total_energy,
                        f"{adsorbate['output']['dir_name']}",
                    )
                else:
                    print(
                        ads_name,
                        total_energy,
                        f"{adsorbate['output']['dir_name']}",
                    )
                #docs['ADS'].append({ads + "*": adsorbate})
                docs[ads+'*'].append(adsorbate)
        if comp_adsorbate == comp_substrate:
            docs['BARE'].append(adsorbate)
            total_energy = f'{adsorbate["output"]["output"]["energy"]:8.2f}'
            print(
                'BARE',
                total_energy,
                f"{adsorbate['output']['dir_name']}",
            )
    return docs


def save_poscar(struct, adsorbate_only=True):
    for name in struct:
        if adsorbate_only and not "*" in name:
            continue
        struct[name].to_file(f'POSCAR.{name.replace("*","-")}')
