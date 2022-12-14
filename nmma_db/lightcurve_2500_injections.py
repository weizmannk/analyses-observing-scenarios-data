import numpy as np
import os
import pandas as pd
from astropy.time import Time
from os.path import exists
from pathlib import Path
from tqdm.auto import tqdm
from astropy.time import Time
import time
from datetime import timedelta

#this directory need to be absolute dir (star from user home)
outdir = '/home/weizmann.kiendrebeogo/NMMA/nmma/observing-scenarios-lc/lc-2500-injection'
output ='./lightcurve'
if not os.path.isdir(output):
    os.makedirs(output)

#binary outdir type 
#dir_type  = ['outdir_BNS', 'outdir_NSBH']

# telescope names 
telescopes = ['ztf', 'rubin']
run_names = ['O4', 'O5']
pops  = ['NSBH']

seed_list=[816, 323, 364, 564, 851]
#exposure only ztf case
exposure = [180, 300]
#number total of injection
n_inj = 2500


print("Start : %s" % time.ctime())
start_time = time.time()
time.sleep(3)

with tqdm(total=len(telescopes) * len(run_names)*len(pops)) as progress:
    for telescope in tqdm(telescopes):
        for run_name in run_names:
            for pop in pops:
                job, s = 0, 0
                path = Path(f'{outdir}/{telescope}/{run_name}/outdir_{pop}/{pop}')
                lcs=pd.DataFrame()
                if telescope == 'ztf':

                    for ii in range(1, int(n_inj)+1):
                        for exp in exposure:
                            for seed in seed_list:
                                sim= Path(f'{path}/{str(ii)}_{str(exp)}_{str(seed)}')
                                try:
                                    if exists(f'{sim}/lightcurves.png'):
                                        print('yes' )
                                        lc=pd.read_csv(f'{sim}/lc.csv')
                                        lc.sort_values(by='jd', ascending=True)
                                        lc['sim']=sim
                                        lcs=pd.concat([lcs,lc], ignore_index = True)
                                        print('success on : ', sim)
                                        s+=1
                                except FileNotFoundError as e:
                                    print(f'{sim} missing')
                                job+=1

                else:
                    if telescope == 'rubin':
                        for ii in range(1, int(n_inj)+1):
                            for seed in seed_list:
                                sim= Path(f'{path}/{str(ii)}_{str(seed)}')
                                try:
                                    if exists(f'{sim}/lightcurves.png'):
                                        lc=pd.read_csv(f'{sim}/lc.csv')
                                        lc.sort_values(by='jd', ascending=True)
                                        lc['sim']=sim
                                        lcs=pd.concat([lcs,lc], ignore_index = True)
                                        print('success on : ', sim)
                                        s+=1
                                except FileNotFoundError as e:
                                    print(f'{sim} missing')
                                job+=1

                #time is in mjd and the time column name in jd,
                #so to ovoid any confusion we convert the time values in jd 
                lcs['jd'] = Time(lcs['jd'], format='mjd').jd

                nmma_dict = { 'jd':'jd', 'filter': 'filter', 'mag':'mag',  'mag_unc': 'mag_unc',   'limmag': 'limmag', 'programid': 'programid', 'sim': 'sim'}

                lcs.rename(columns=nmma_dict, inplace=True)
                lcs.reset_index(drop=True)
                lcs.pop('Unnamed: 0')
                lcs.sort_values(by='jd', ascending=True)
                lcs.to_csv(f'{output}/{telescope}_{run_name}_{pop}.csv', index=False)
                del lcs, nmma_dict
                progress.update()


                print("=====================================================")
                print("    ")
                print("Total number of jobs = ", job)
                print('The number of job where is at least one detection :', s)
                print("    ")
                print("=====================================================")
                print("    ")
                print("    ")
                time.sleep(3)

print("End : %s" % time.ctime())
time.sleep(2)

# determine the time take to run this, in the format (Hours,Minutes, and Secondes)
hours, rem       =  divmod((time.time() - start_time), 3600)
minutes, seconds = divmod(rem, 60)


print("    ")
print("    ")

print("+----------------------------------------------------+")
print("    ")
print("--- This run took {:0>2}:{:0>2}:{:05.2f} ---".format(int(hours),int(minutes),seconds))
print("    ")     
print("+----------------------------------------------------+")
