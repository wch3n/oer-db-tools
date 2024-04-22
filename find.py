#!/usr/bin/env python3

''' queries the mongo db for a series of OER intermediates on a specific substrate;
    results are stored in the docs dictionary for the gas phase molecules, surface, 
    and all intermediates found from the db'''

import oer

store = oer.connect_db()
substrate = 'Mo4 C3 O2'
substrate = 'Mo Ti Nb V C3 O2' 
#substrate = 'Mo12 Ti12 Nb12 V12 C36 O23' 
functional = 'PBE'
docs = oer.find_all(store, substrate, functional, reaction='HER', series='/her/')
