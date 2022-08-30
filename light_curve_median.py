import os
import numpy as np
from pathlib import Path
from astropy.io import ascii
import pandas as pd
from tqdm.auto import tqdm
import matplotlib.pyplot as plt

path_dir =  f"{os.path.dirname(os.path.realpath('__file__'))}/lightcurve"

outdir = "./outdir"
if not os.path.isdir(outdir):
    os.makedirs(outdir)

run_names  = ['O4'] #, 'O5']
distribution = ['Petrov'] #, 'Farah']
telescopes  = ['ztf', 'rubin'] 
pops       = ['BNS', 'NSBH']

with tqdm(total=len(distribution)*len(telescopes) * len(run_names)*len(pops)) as progress:
    for run_name in run_names:  
        for dist in distribution:
            plt.clf()
            # Figure Plot
            fig, axs = plt.subplots(nrows=2, ncols=2) 
            for telescope in telescopes:   
                if telescope == "zft":
                    ax1, ax2 = 0, 0 
                else:
                    ax1, ax2 = 0, 1 
                            
                for pop in pops:     
                    path = Path(f"{path_dir}/{dist}_{telescope}_{run_name}_{pop}.csv")
                    lc = ascii.read(f"{path}", format='csv')
                    lc = lc[lc['mag_unc']!=99.0] 
                    
                    for filt in ["g", "r", "i"] : 
                        tab = lc[(lc['filter'] == filt)]
                        tab.sort('jd')
                        mag = tab['mag']
                        
                        if filt == 'g':
                            c ='g'
                        elif filt == 'r':
                            c ='r'
                        elif filt == 'i':
                            c = 'b'
                        else:
                            print("This filter is no take care in this analysis")
                        
                        axs[ax1, ax2].hist(mag, bins=200, density =12, histtype='step' , color=c,  label= filt , linewidth=1.3)

                    axs[ax1, ax2].legend(prop={'size': 10})
                    axs[ax1, ax2].set_title(f'{dist} {telescope} {run_name}', fontname="Times New Roman", size=13, fontweight="bold") 
                    ax1, ax2 = 1, ax2
                    
                    progress.update()
        
            plt.gcf().set_size_inches(12, 12)
            plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.4,
                    hspace=0.4)
            fig.tight_layout()
            plt.savefig(f'{outdir}/{dist}_magnitude_{run_name}.png')
            plt.close()

