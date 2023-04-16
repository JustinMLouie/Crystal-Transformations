import json
import numpy as np
from pymatgen.core import Structure
from sys import argv

data = [json.loads(line) for line in open('db.json.gz', 'r')]

for i in range(len(data)):
    entry = data[i]

    if entry['material_id'] == argv[1]:
        struc = Structure.from_dict(entry['structure'])
        struc.to(filename='POSCAR_{}'.format(argv[2]))
