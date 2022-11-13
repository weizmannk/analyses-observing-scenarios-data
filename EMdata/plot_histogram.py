import os
import numpy as np
from pathlib import Path
from astropy.io import ascii
import pandas as pd
from tqdm.auto import tqdm
import matplotlib.pyplot as plt

path_dir =  f"{os.path.dirname(os.path.realpath('__file__'))}/output_lc"

outdir = "./output_lc"
if not os.path.isdir(outdir):
    os.makedirs(outdir)

run_names  = ['O4']
distribution = ['Farah']
telescopes  = ['rubin'] 
pops       = ['NSBH']

median_tables = {}
with tqdm(total=len(distribution)*len(telescopes) * len(run_names)*len(pops)) as progress:
    for dist in distribution:
        median_tables[dist] = {}
        for run_name in run_names:
            median_tables[dist][run_name] = {}
            plt.clf()
            # Figure Plot
            # remove the number of columns
            fig, axs = plt.subplots(nrows=1, ncols=1) 
          
            for telescope in telescopes:   
                median_tables[dist][run_name][telescope] = {}
                for pop in pops:   
                    median_tables[dist][run_name][telescope][pop] = {}  
                    path = Path(f"{path_dir}/{dist}_lc_{telescope}_{run_name}_{pop}.csv")
                    lc = ascii.read(f"{path}", format='csv')
                    lc = lc[lc['mag_unc']!=99.0]  
                    
                    for filt in ["g", "r", "i"] : 
                        tab = lc[(lc['filter'] == filt)]
                        tab.sort('jd')
                        mag = tab['mag']
                        mag_err = tab['mag_unc']
                        
                        if len(tab) !=0:
                            median_tables[dist][run_name][telescope][pop][filt] = f'{np.median(mag)} +/- {np.median(mag_err)}'
                                
                            if filt == 'g':
                                c ='g'
                            elif filt == 'r':
                                c ='r'
                            elif filt == 'i':
                                c = 'k'
                                
                            else:
                                print("This filter is no take care in this analysis")
                             
                            if telescope == 'rubin':
                                if pop == 'NSBH': 
                                    axs.hist(mag, bins=5, density =12, histtype='step' , color=c,  label= filt , linewidth=1.2)
                                    axs.legend(prop={'size': 10}, loc='lower left')
                                    #axs.set_title(f'{run_name} {dist} {telescope} ', fontname="Times New Roman", size=13, fontweight="bold")
                                    axs.text(0.03, 0.95, f'{run_name} {dist}  {telescope} ', transform=axs.transAxes,  va='top', fontname="Times New Roman", size=13, fontweight="bold")
                                    axs.text(0.05, 0.85, pop, transform=axs.transAxes, color='blue',  va='top', fontname="Times New Roman", size=13, fontweight="bold")
                                    axs.set_xlabel(r'Magnitude', fontname="Times New Roman", size=13, fontweight="bold")  
                    

                    progress.update()
            
            plt.gcf().set_size_inches(6, 6)
            plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.4,
                    hspace=0.4)
            fig.tight_layout()
            plt.savefig(f'{outdir}/{dist}_magnitude_{run_name}_{pop}.png')
            plt.close()
