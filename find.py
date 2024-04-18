#!/usr/bin/env python3

''' queries the mongo db for a series of OER intermediates on a specific substrate;
    results are stored in the docs dictionary for the gas phase molecules, surface, 
    and all intermediates found from the db'''

import oer

store = oer.connect_db()
substrate = 'Ti4 C3 O2'
#substrate = 'Mo Ti Nb V C3 O2' 
functional = 'PBE'
docs = oer.find_all(store, substrate, functional)

# slab for example
doc = docs['SLAB']['output']['output']
print(doc.keys())
print(doc['energy'])

# OH* adsorbate
doc = docs['OH*']['output']['output']
print(doc['energy'])
