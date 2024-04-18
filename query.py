#!/usr/bin/env python3

import oer

store = oer.connect_db()
#substrate = 'Ti4 C3 O2'
substrate = 'Mo Ti Nb V C3 O2'
functional = 'PBE'
anchor_index = 31
deltaG, energy, struct = oer.report(store, substrate, functional, anchor_index=anchor_index, write_poscar=True)
print(deltaG, energy)
