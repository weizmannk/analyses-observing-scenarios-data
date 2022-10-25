import os
import numpy as np
from statistics import median, mean
from pathlib import Path
from astropy.io import ascii
import pandas as pd
from tqdm.auto import tqdm

path_dir =  f"{os.path.dirname(os.path.realpath('__file__'))}/lightcurve"

outdir = "./outdir"
if not os.path.isdir(outdir):
    os.makedirs(outdir)

run_names  = ['O4', 'O5']
distribution = ['Petrov_lc', 'Farah_lcs']
telescopes  = ['ztf', 'rubin'] 
pops       = ['BNS', 'NSBH']

tables = {}
with tqdm(total=len(distribution)*len(telescopes) * len(run_names)*len(pops)) as progress:
    for run_name in run_names: 
        tables[run_name] = {} 
        for dist in distribution:
            tables[run_name][dist] = {}
            for telescope in telescopes: 
                tables[run_name][dist][telescope] = {}    
                for pop in pops:
                    tables[run_name][dist][telescope][pop] = {}    
                    path = Path(f"{path_dir}/{dist}_{telescope}_{run_name}_{pop}.csv")
                    lc = ascii.read(f"{path}", format='csv')
                    lc = lc[lc['mag_unc']!=99.0] 
                    
                    for filt in ["g", "r", "i"] : 
                        tab = lc[(lc['filter'] == filt)]
                        tab.sort('jd')
                        mag = tab['mag']
                        mag_err = tab['mag_unc']
                        
                        tables[run_name][dist][telescope][pop][filt] = median(mag)
                        tables[run_name][dist][telescope][pop][f'{filt}_err'] = mean(mag_err)

                    progress.update() 
