import os
from astropy.table import Table
import pandas as pd

infile = f"{os.path.dirname(os.path.realpath(__file__))}/O1O2O3all_mass_h_iid_mag_iid_tilt_powerlaw_redshift_maxP_events_all.h5"

data = Table.read(infile)

outfile ='farah.h5'
Table({
    'mass1': data['mass_1'],
    'mass2': data['mass_2']
    }).write(outfile, overwrite=True)
 

df = pd.DataFrame({
    'mass1': data['mass_1'],
    'mass2': data['mass_2']
    }
).sort_values("mass1").reset_index(drop=True)


if not os.path.isfile('farah.tex'):
    with open('farah.tex', 'w') as f:
        f.write(df.to_latex(
            longtable=True, 
            caption='sThe farah population (BNS, NSBH, BBH)', 
            label='tab:farah pop')
        )

else:
    os.remove('farah.tex')
    with open('farah.tex', 'w') as f:
        f.write(df.to_latex(
            longtable=True, 
            caption='The farah population (BNS, NSBH, BBH)', 
            label='tab:fararah pop')
        )



print(f'Farah Method :')
print('=====================================================')
print(f"number of BNS: {sum(df['mass1']<3)}, number of NSBH:"
f"{sum(df['mass2']<3) - sum(df['mass1']<3)}, number of BBH: "
f"{sum(df['mass2']>3)}")

print('      \n')

print('=====================================================')

print('Weizmann Process')


print('=====================================================')
print(f"number of BNS: {len(df[(df['mass1']<3)])}, number of NSBH:"
f"{len(df[(df['mass1']>=3) & (df['mass2']<3)])}, number of BBH: "
f"{len(df[(df['mass2']>=3)])}")

print('=====================================================')
