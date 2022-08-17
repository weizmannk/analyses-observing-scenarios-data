import os
from astropy.table import Table
import pandas as pd

infile = f"{os.path.dirname(os.path.realpath(__file__))}/O1O2O3all_mass_h_iid_mag_iid_tilt_powerlaw_redshift_maxP_events_all.h5"

outfile ='farah.tex' 
data = Table.read(infile)

df = pd.DataFrame({
    'mass1': data['mass_1'],
    'mass2': data['mass_2'],
        }
    ).sort_values("mass1")	


 
if not os.path.isfile(outfile):
    with open('farah.tex', 'w') as f:
        f.write(df.to_latex(
            longtable=True, 
            caption='sThe farah population (BNS, NSBH, BBH)', 
            label='tab:farah pop')
        )

else:
    os.remove(outfile)
    with open('farah.tex', 'w') as f:
        f.write(df.to_latex(
            longtable=True, 
            caption='The farah population (BNS, NSBH, BBH)', 
            label='tab:fararah pop')
        )

print(
    df.to_latex(
        longtable=True, 
        caption='The farah population (BNS, NSBH, BBH)', 
        label='tab:fararah pop'
        )
    )


print( "=====================================")

print(len(df))

print( "=====================================")



"""
with open("farah.tex", "w") as f:
    f.write("\\begin{tabular}{" + " | ".join(["c"] * len(df.columns)) + "}\n")
    for i, row in df.iterrows():
        f.write(" & ".join([str(x) for x in row.values]) + " \\\\\n")
    f.write("\\end{tabular}")
"""









