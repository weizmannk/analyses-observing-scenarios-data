import os
import numpy as np
from pathlib import Path
from astropy.io import ascii
import pandas as pd
from tqdm.auto import tqdm
import matplotlib.pyplot as plt
#plt.style.use('seaborn-v0_8-paper')
import pandas as pd

import random

import matplotlib as mpl
mpl.use('agg')
import matplotlib.gridspec as gridspec




import seaborn as sns 



#sns.set_style("darkgrid")


path_dir =  f"{os.path.dirname(os.path.realpath('__file__'))}/lightcurve/output_lc"

outdir = "./paper_plots"
if not os.path.isdir(outdir):
    os.makedirs(outdir)

ztf_filters = ["g", "r", "i"]
rubin_filters = ["g","r","i","u","z","y"]

run_names  = ['O4']
distribution = ['Farah']
telescopes  = ['ZTF', 'Rubin'] 
pops       = ['BNS']

with tqdm(total=len(distribution)*len(telescopes) * len(run_names)*len(pops)) as progress:
    for dist in distribution:
        for run_name in run_names:
            plt.clf()
            # Figure Plot
            fig, axs = plt.subplots(nrows=2, ncols=2)
            
                    #fig = plt.figure(figsize=(6,6))
            #fig.subplots_adjust(left=0.02,right=0.98,bottom=0.1,top=0.98,wspace=0.3,hspace=0.3)
            rows,cols=2,2
            gs = gridspec.GridSpec(rows,cols)
            sax = []
            
            for r in range(rows):
                for c in range(cols):
                    sax.append(plt.subplot(gs[cols*r+c]))
            
            for telescope in telescopes:   
                #median_tables[dist][run_name][telescope] = {}
                for pop in pops:   
                    #median_tables[dist][run_name][telescope][pop] = {}  
                    path = Path(f"{path_dir}/{dist}_lc_{telescope}_{run_name}_{pop}.csv")
                    #lc = ascii.read(f"{path}", format='csv')
                    lc = pd.read_csv(f"{path}")
                    lc=lc[lc['mag_unc']!=99.0]
                    #split_inj =np.unique(lc['sim'])
                    
                    #for inj in split_inj:
                     #   indx = np.where(lc['sim'] == inj)
                        
                      #  lc = lc[lc[indx]]

                    lc.groupby('sim').ngroups

                    band_colors={'g':'limegreen','r':'orangered','i':'goldenrod', 'u': 'blue', 'y':'darkviolet', 'z':'k' }#,'limegreen','darkturquoise',
                    

                    
                    hist, bins = np.histogram(lc.groupby(['sim']).apply(lambda x: x.shape[0]),bins=np.arange(-0.5,10.,1))
                        
                    if telescope =='ZTF':

                        for name,group in lc.groupby(['sim','filter']).apply(lambda x: x['mag'].min()).groupby(level=1):
                            hist, bins = np.histogram(group,bins=np.arange(16,22,0.5))
                            sax[0].step(bins, np.pad(hist, (1, 0)) / hist.max(), alpha=0.5,color=band_colors[name],lw=2)
                        hist, bins = np.histogram(lc.groupby(['sim']).apply(lambda x: x['mag'].min()),bins=np.arange(16,22,0.5))
                        sax[0].step(bins, np.pad(hist, (1, 0)) / hist.max(),color='k',lw=0.5)

                        sax[0].set_xlabel('peak mag')
                        sax[0].get_yaxis().set_visible(False)

                        sax[1].scatter(lc.groupby(['sim']).apply(lambda x: x.shape[0]),lc.groupby(['sim']).apply(lambda x: x['mag'].min()),s=1)
                        sax[1].set_xlabel('num photo (all)')
                        sax[1].set_ylabel('peak mag (all)')
                        sax[1].invert_yaxis()
                    
                    else:

                        for name,group in lc.groupby(['sim','filter']).apply(lambda x: x['mag'].min()).groupby(level=1):
                            hist, bins = np.histogram(group,bins=np.arange(16,22,0.5))
                            sax[2].step(bins, np.pad(hist, (1, 0)) / hist.max(), alpha=0.5,color=band_colors[name],lw=2)
                        hist, bins = np.histogram(lc.groupby(['sim']).apply(lambda x: x['mag'].min()),bins=np.arange(16,22,0.5))
                        sax[2].step(bins, np.pad(hist, (1, 0)) / hist.max(),color='k',lw=0.5)

                        sax[2].set_xlabel('peak mag')
                        sax[2].get_yaxis().set_visible(False)

                        sax[3].scatter(lc.groupby(['sim']).apply(lambda x: x.shape[0]),lc.groupby(['sim']).apply(lambda x: x['mag'].min()),s=1)
                        sax[3].set_xlabel('num photo (all)')
                        sax[3].set_ylabel('peak mag (all)')
                        sax[3].invert_yaxis()
                        
                        progress.update()
            
        plt.gcf().set_size_inches(9, 6)
        plt.subplots_adjust(left=0.1,
                bottom=0.1,
                right=0.9,
                top=0.9,
                wspace=0.4,
                hspace=0.4)
        fig.tight_layout()
        plt.savefig(f'{dist}_magnitude_{run_name}_{telescope}_{pop}.png',dpi=300)
        plt.close()

