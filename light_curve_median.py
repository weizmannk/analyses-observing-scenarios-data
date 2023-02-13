import os
import numpy as np
from pathlib import Path
from astropy.io import ascii
import pandas as pd
from tqdm.auto import tqdm
import matplotlib.pyplot as plt
#plt.style.use('seaborn-v0_8-paper')
import seaborn as sns 

sns.set_style("darkgrid")


path_dir =  f"{os.path.dirname(os.path.realpath('__file__'))}/lightcurve/output_lc"

outdir = "./paper_plots"
if not os.path.isdir(outdir):
    os.makedirs(outdir)

ztf_filters = ["g", "r", "i"]
rubin_filters = ["g","r","i","u","z","y"]

run_names  = ['O4', 'O5']
distribution = ['Farah']
telescopes  = ['ztf', 'rubin'] 
pops       = ['BNS', 'NSBH']

median_tables = {}
with tqdm(total=len(distribution)*len(telescopes) * len(run_names)*len(pops)) as progress:
    for dist in distribution:
        
        median_tables[dist] = {}
        for run_name in run_names:
            median_tables[dist][run_name] = {}
            plt.clf()
            # Figure Plot
            fig, axs = plt.subplots(nrows=2, ncols=2) 
          
            for telescope in telescopes:   
                median_tables[dist][run_name][telescope] = {}
                for pop in pops:   
                    median_tables[dist][run_name][telescope][pop] = {}  
                    path = Path(f"{path_dir}/{dist}_lc_{telescope}_{run_name}_{pop}.csv")
                    lc = ascii.read(f"{path}", format='csv')
                    lc = lc[lc['mag_unc']!=99.0]
                    
                    if telescope== 'ztf':
                        filters = ztf_filters
                    else :
                        filters = rubin_filters                          

                    for filt in filters : 
                        tab = lc[(lc['filter'] == filt)]
                        tab.sort('jd')
                        mag = tab['mag']
                        mag_err = tab['mag_unc']
                        
                        if len(tab) !=0:
                            median_tables[dist][run_name][telescope][pop][filt] = f'{np.median(mag)} +/- {np.median(mag_err)}'
                                                     
                            if filt == "g":
                                c ="darkgreen"
                            elif filt == "r":
                                c ="darkred"
                            elif filt == "i":
                                c = "orange"
                            elif filt == "u":
                                c ="blue"
                            elif filt == "z":
                                c ="k"
                            elif filt == "y":
                                c = "darkviolet" 
                            
                            else:    
                                print("This filter is no take care in this analysis")
                          
                            if telescope == 'ztf':
                                if pop == 'BNS': 
                                    axs[0, 0].hist(mag, bins=70, histtype='step' , color=[c],  label= filt , linewidth=1.2)
                                    axs[0, 0].legend(prop={'size':10}) 
                                    #axs[0, 0].set_title(f'{run_name} {dist}  {telescope} ', fontname="Times New Roman", size=14, fontweight="bold")
                                
                                    axs[0, 0].text(0.03, 0.95, r'ZTF', transform=axs[0,0].transAxes,  va='top', fontname="Times New Roman", size=14, fontweight="bold")
                                    axs[0, 0].text(0.03, 0.85, pop, transform=axs[0,0].transAxes, color='green',  va='top', fontname="Times New Roman", size=14, fontweight="bold")
                                    axs[0, 0].set_xlabel(r'Magnitude', fontname="Times New Roman", size=14, fontweight="bold") 
                                    
                                else:
                                    axs[1, 0].hist(mag, bins=50, histtype='step' , color=[c],  label= filt , linewidth=1.2)
                                    #axs[1, 0].set_title(f'{run_name} {dist} {telescope} ', fontname="Times New Roman", size=14, fontweight="bold")
                                    axs[1, 0].text(0.03, 0.95, r'ZTF', transform=axs[1,0].transAxes,  va='top', fontname="Times New Roman", size=14, fontweight="bold")
                                    axs[1, 0].text(0.03, 0.85, pop, transform=axs[1,0].transAxes, color='blue',  va='top', fontname="Times New Roman", size=14, fontweight="bold")
                                    axs[1, 0].set_xlabel(r'Magnitude', fontname="Times New Roman", size=14, fontweight="bold")  
                                
                            else:
                                if pop == 'BNS': 
                                    axs[0, 1].hist(mag, bins=50, histtype='step' , color=[c],  label= filt , linewidth=1.2)
                                    #axs[0, 1].set_title(f'{run_name} {dist} {telescope} ', fontname="Times New Roman", size=14, fontweight="bold")
                                    axs[0, 1].text(0.03, 0.95, r'LSST', transform=axs[0,1].transAxes,  va='top', fontname="Times New Roman", size=14, fontweight="bold")
                                    axs[0, 1].text(0.03, 0.85, pop, transform=axs[0,1].transAxes, color='green',  va='top', fontname="Times New Roman", size=14, fontweight="bold")
                                    axs[0, 1].set_xlabel(r'Magnitude', fontname="Times New Roman", size=14, fontweight="bold") 
                                else:
                                    axs[1, 1].hist(mag, bins=25, histtype='step' , color=[c],  label= filt , linewidth=1.2)
                                    axs[1, 1].legend(prop={'size':10}, loc='lower left')
                                    #axs[1, 1].set_title(f'{run_name} {dist} {telescope} ', fontname="Times New Roman", size=14, fontweight="bold")
                                    axs[1, 1].text(0.03, 0.95, r'LSST', transform=axs[1,1].transAxes,  va='top', fontname="Times New Roman", size=14, fontweight="bold")
                                    axs[1, 1].text(0.03, 0.85, pop, transform=axs[1,1].transAxes, color='blue',  va='top', fontname="Times New Roman", size=14, fontweight="bold")
                                    axs[1, 1].set_xlabel(r'Magnitude', fontname="Times New Roman", size=14, fontweight="bold")  
                    
                    progress.update()
            
            plt.gcf().set_size_inches(9, 6)
            plt.subplots_adjust(left=0.1,
                    bottom=0.1,
                    right=0.9,
                    top=0.9,
                    wspace=0.4,
                    hspace=0.4)
            fig.tight_layout()
            plt.savefig(f'{outdir}/{dist}_magnitude_{run_name}.png')
            plt.close()
