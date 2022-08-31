from scipy import stats
from astropy.table import Table, vstack
import numpy as np

n_samples = 1000000

# Atrophysique Binaries Black Hole  
# minimal and  max spin
bh_astro_spin_min = -0.99
bh_astro_spin_max = +0.99
ns_astro_spin_min = -0.05
ns_astro_spin_max = +0.05

ns_mass_min = 1.0
ns_mass_max = 3.0
bh_mass_min = 3.0
bh_mass_max = 50.0

# ns mass_dist
ns_astro_mass_distn = stats.norm(1.33, 0.09)
ns_astro_spin_distn = stats.uniform(ns_astro_spin_min, ns_astro_spin_max - ns_astro_spin_min)

# bh mass_distn
bh_astro_mass_distn = stats.pareto(b=1.3)
bh_astro_spin_distn = stats.uniform(bh_astro_spin_min, bh_astro_spin_max - bh_astro_spin_min)

# censor the mass distribution, remove or ovoid a non BNS-astro, NSBH-astro and  BBH-astro mass 
def draw_masses(n_samples, mass_distn, mass_min, mass_max):
    nbad = n_samples
    mass = np.empty(n_samples)
    bad = np.ones(n_samples, dtype=bool)
    while nbad > 0:
        mass[bad] = mass_distn.rvs(nbad)
        bad = (mass < mass_min) | (mass > mass_max)
        nbad = np.sum(bad)
    return mass

pops = [ "BNS", "BBH", "NSBH"]
tables = {}

for pop in pops:
    
    tables[pop] = {}

    if pop =="BNS":
        mass1  = draw_masses(n_samples, ns_astro_mass_distn, ns_mass_min, ns_mass_max)
        mass2  = draw_masses(n_samples, ns_astro_mass_distn, ns_mass_min, ns_mass_max)
        spin1z = ns_astro_spin_distn.rvs(n_samples)
        spin2z = ns_astro_spin_distn.rvs(n_samples)
        
        
    elif pop =="BBH":
        mass1  = draw_masses(n_samples, bh_astro_mass_distn, bh_mass_min, bh_mass_max)
        mass2  = draw_masses(n_samples, bh_astro_mass_distn, bh_mass_min, bh_mass_max)
        spin1z = bh_astro_spin_distn.rvs(n_samples)
        spin2z = bh_astro_spin_distn.rvs(n_samples)
    
   # NSBH 
    else:
        mass1  = draw_masses(n_samples, bh_astro_mass_distn, bh_mass_min, bh_mass_max)
        mass2  = draw_masses(n_samples, ns_astro_mass_distn, ns_mass_min, ns_mass_max)
        spin1z = bh_astro_spin_distn.rvs(n_samples)
        spin2z = ns_astro_spin_distn.rvs(n_samples)
  

    # swap masses to ensure that mass1 >= mass2 
    swap = mass1 < mass2
    mass1[swap], mass2[swap] = mass2[swap].copy(), mass1[swap].copy()
    # We could simply use this one swap
    #mass1, mass2 = np.maximum(mass1, mass2), np.minimum(mass1, mass2)


    tables[pop]['mass1'] = mass1
    tables[pop]['mass2'] = mass2
    tables[pop]['spin1z'] = spin1z
    tables[pop]['spin2z']  = spin2z



bns_astro   = Table(tables['BNS'])
nsbh_astro  = Table(tables['NSBH'])
bbh_astro   = Table(tables['BBH'])


try:
    Petrov =   vstack([bns_astro, nsbh_astro, bbh_astro], join_type='exact')
        
except TableMergeError as ex:
    print(ex)

else:
    # save data on .h5 file 
    Petrov.write(
    f"./data/Pretrov_Farah_input/petrov.h5", overwrite=True
)

 

data = Table.read("./data/Pretrov_Farah_input/petrov.h5")
