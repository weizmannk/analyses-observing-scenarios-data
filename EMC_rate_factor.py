import numpy as np 
from astropy.table import join, Table
from scipy import stats



# Lower 5% and upper 95% quantiles of log normal distribution
rates_table = Table(
    [
        # O3 R&P paper Table II row 1 last column
        {'population': 'BNS', 'lower': 100., 'mid': 240., 'upper': 510.},
        {'population': 'NSBH', 'lower': 100., 'mid': 240., 'upper': 510.},
        {'population': 'BBH', 'lower': 100., 'mid': 240., 'upper': 510.}
    ]
)

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

standard_90pct_interval, = np.diff(stats.norm.interval(0.9))
rates_table['mu'] = np.log(rates_table['mid'])
rates_table['sigma'] = (np.log(rates_table['upper']) - np.log(rates_table['lower'])) / standard_90pct_interval

rates_table

#fiducial_log_rates = np.asarray(rates_table['mu'])
#fiducial_log_rate_errs = np.asarray(rates_table['sigma'])


