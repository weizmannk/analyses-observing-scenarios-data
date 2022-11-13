1################################

To split NSBH from injection file O4 
First osf of clone this reposetory in 'EMdata'

    ** git clone https://github.com/weizmannk/obs-scenarios-data-2022.git

then please run :

    ** python split_farah_populations.py

Then you will get in this directory NSBH 
injection file : 'farah_post_split/runs/O4/nsbh_farah/injections.dat'

2################################
Use nmma to create NSBH 'Bu2019nsbh_injection.sjon' file:

    ** nmma_create_injection --prior-file ~/NMMA/nmma/priors/Bu2019nsbh.prior --injection-file ./farah_post_split/runs/O4/nsbh_farah/injections.dat --eos-file  ~/NMMA/nmma/example_files/eos/ALF2.dat --binary-type NSBH --n-injection 2500 --original-parameters --extension json --aligned-spin -f ./Farah_lc/outdir_NSBH/Bu2019nsbh_injection


3###############################
Lightcurve analysis using Rubin telescope:

    ** light_curve_analysis_condor --model Bu2019nsbh --prior ~/NMMA/nmma/priors/Bu2019nsbh.prior --svd-path  ~/NMMA/nmma/svdmodels --outdir ./Farah_lc/outdir_NSBH --label injection_Bu2019nsbh --injection ./Farah_lc/outdir_BNS/Bu2019lnsbh_injection.json --injection-num 2500 --generation-seed 816 323 364 564 851  --condor-dag-file condor.dag --condor-sub-file condor.sub --bash-file condor.sh

### the jobs will be save here './Farah_lc/outdir_NSBH/NSBH'

4##############################
Submit jobs

** condor_submit_dag -f condor.dag


5############################
Extract light curve csv data from each 
job where there is at least one detection

thenrun this file 

    ** python light_curve.py

the you will get an csv file in this directory 'output_lc/Farah_lc_rubin_O4_NSBH.csv'
where is the all detection together.


6################################
plot histogramme for the light curve magnitude

then run this file:

    ** python plot_histogram.py

you will get in this directory 'output_lc/Farah_magnitude_O4_NSBH.png'