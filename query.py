#!/usr/bin/env python3

import oer
import pandas as pd

store = oer.connect_db()
substrate = 'V4 C3 O2'
#substrate = 'Mo Ti Nb V C3 O2'
substrate = 'Mo12 Ti12 Nb12 V12 C36 O23' 
functional = 'PBE'

anchors = [6,7,8,18,19,20,30,31,32,42,43,44]
anchors = [[6,7,8]]
name = 'o93'
out = {}
for anchor in anchors:
    deltaG, energy, struct = oer.report(store, substrate, functional, anchor_index=anchor, write_poscar=False, series=name)
    out[name] = deltaG

print(out)
df = pd.DataFrame.from_dict(out)
df.to_csv('{0}_{1}.csv'.format(substrate.lower().replace(' ',''), functional.lower()))
