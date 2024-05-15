#!/usr/bin/env python3

import oer

store = oer.connect_db()

# the quaternary high-entropy mxene
# for pure mxenes, substrate = "Ti4 C3 O2", "Mo4 C3 O2", ...
substrate = 'Mo12 Ti12 Nb12 V12 C36 O24'

# PBE (no dispersion)
functional = 'PBE'

# list of adsorption sites for OER on the high-entropy mxene
# for each element there are three surface sites
sites_Ti = [6,7,8]
sites_Nb = [18,19,20]
sites_V = [30,31,32]
sites_Mo = [42,43,44]

# for pure mxenes, the adsorption site index is always 24

# let's check the Ti adsorption site with index 6:
for anchor in sites_Ti[:1]:
    deltaG, energy, struct, forces = oer.report(store, substrate, functional, anchor_index=anchor, write_poscar=False, series=None)

# about the results:
# deltaG is a dict of the Gibbs free energies along the OER reaction coordinates (1 to 4).
# energy is a dict where the total energies of the various adsorbates are stored.
print(energy.keys())

# struct and forces share the same data strucure as energy
# to get structure of the substrate 
print(struct['substrate'])

# the structure can be passed to ASE via `AseAtomsAdaptor` of pymatgen

# to get forces of the OH* adsorbate on the substrate 
print(forces['OH*'])
