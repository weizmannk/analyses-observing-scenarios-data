import pandas as pd
import matplotlib.pyplot as plt
import sys
import argparse
import json

#parser=argparse.ArgumentParser()
#parser.add_argument("typ", type=str)
#parser.add_argument("n_inj", type=int)
#parser.add_argument("--seed_list", nargs="+", type=int)

outdir='outdir'
#args=parser.parse_args()

seed_list=[816, 323, 364, 564, 851]
n_inj =2500
type  = ['BNS']

for bin_type in type:
   
  with open(f"{outdir}/{bin_type}/injection.json", "r") as read_content:
    inj=json.load(read_content)

  lcs=pd.DataFrame()
  for i in range(n_inj):
    for exp in [180,300]:
      for seed in seed_list:
        sim=f"{outdir}{bin_type}/{i}_{exp}_{seed}"
        dL=inj['injections']['content']['luminosity_distance'][i]
        try:
          lc=pd.read_csv(sim+'/too.csv')
          lc['sim']=sim
          lc['luminosity_distance']=dL
          lcs=pd.concat([lcs,lc])  
        except FileNotFoundError as e:
          print(sim+' missing')

  lcs=lcs.reset_index(drop=True)
  lcs.to_csv(f"{bin_type}.csv",index=False)

  lcs=pd.DataFrame()
  for i in range(n_inj):
    for seed in seed_list:
      sim=f"{outdir}{bin_type}/{i}_{exp}_{seed}_all"
      dL=inj['injections']['content']['luminosity_distance'][i]
      try:
        lc=pd.read_csv(sim+'/lc.csv')
        lc['sim']=sim
        lc['luminosity_distance']=dL
        lcs=pd.concat([lcs,lc[['jd','mag','filter','sim','luminosity_distance']]])
      except FileNotFoundError as e:
        print(sim+' missing')

  lcs=lcs.reset_index(drop=True)
  lcs.to_csv(f"{bin_type}_all.csv",index=False)

