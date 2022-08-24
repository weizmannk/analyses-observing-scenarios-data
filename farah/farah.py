#!/usr/bin/env python
"""Convert Farah's distribution to suitable format for bayestar-inject."""
import os
from astropy.table import Table

input =f"{os.path.dirname(os.path.realpath(__file__))}/O1O2O3all_mass_h_iid_mag_iid_tilt_powerlaw_redshift_maxP_events_all.h5"

output = os.path.abspath("..") + "/" + os.path.join("data/farah_input")
if not os.path.isdir(output):
    os.makedirs(output)

 
data = Table.read(input)
Table({
    'mass1': data['mass_1'],
    'mass2': data['mass_2'],
    'spin1z': data['a_1'] * data['cos_tilt_1'],
    'spin2z': data['a_2'] * data['cos_tilt_2']
}).write(
    f"{output}/farah.h5",
    overwrite=True
)

# print the number of each population in Farah distribution
print(f"number of BNS: {len(data[(data['mass_1']<3)])}, number of NSBH:"
f"{len(data[(data['mass_1']>=3) & (data['mass_2']<3)])}, number of BBH: "
f"{len(data[(data['mass_2']>=3)])}")

