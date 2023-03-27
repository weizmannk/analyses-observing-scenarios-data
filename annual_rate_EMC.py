import numpy as np



# For splitting into BNS, NSBH, and BBH populations
ns_max_mass = 3

# Calculate effective rate density for each sub-population
table = Table.read('farah.h5')
source_mass1 = table['mass1']
source_mass2 = table['mass2']
rates_table['mass_fraction'] = np.asarray([np.sum((source_mass1 < ns_max_mass) & (source_mass2 < ns_max_mass)),
                                           np.sum((source_mass1 >= ns_max_mass) & (source_mass2 < ns_max_mass)),
                                           np.sum((source_mass1 >= ns_max_mass) & (source_mass2 >= ns_max_mass))]) / len(table)
for key in ['lower', 'mid', 'upper']:
    rates_table[key] *= rates_table['mass_fraction']
del table, source_mass1, source_mass2