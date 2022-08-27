import os
import numpy as np
from pathlib import Path
from astropy.io import ascii
import pandas as pd
from tqdm.auto import tqdm
import matplotlib.pyplot as plt

path_dir =  f"{os.path.dirname(os.path.realpath('__file__'))}/lightcurve"

colors     = ['r', 'g', 'b']

run_names  = ['O4', 'O5']
telescope  = ['ztf', 'rubin'] 
pops       = ['BNS', 'NSBH']

tables = {}

with tqdm(total=len(telescopes) * len(run_names)*len(pops)) as progress:
    for run_name in run_names:
    for telescope in tqdm(telescopes):
       
                    plt.clf()
        # Figure Plot 
        fig, axs = plt.subplots(nrows=3, ncols=2)
            for pop in pops:
                path = Path(f"{path_dir}/{telescope}_{run_name}_{pop}.csv")
                lcs = ascii.read(f"{path}", format='csv')   
                lcs=lcs[lcs['mag_unc']!=99.0]
