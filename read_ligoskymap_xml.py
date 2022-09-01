import os
from tqdm.auto import tqdm
from pathlib import Path

import pandas as pd
from gwpy.table import Table as gwpy_Table

# Read Ligo xml file generate by ligo.skymap
path_dir = "/home/weizmann.kiendrebeogo/OBSERVING_SCENARIOS/observing-scenarios-simulations/runs"

outdir = "./data/inital-distribution-in-bayestar-inject_xml"

if not os.path.isdir(outdir):
    os.makedirs(outdir)
    
def xml_to_dataframe(prior_file):
    table = gwpy_Table.read(prior_file, format="ligolw", tablename="sim_inspiral")
    injection_values = {
        "mass1": [],
        "mass2": [],
        "distance": [],
        "psi": [],
        "phase": [],
        "geocent_time": [],
        "ra": [],
        "dec": [],
        "spin1z": [],
        "spin2z": [],
    }
    for row in table:
        injection_values["mass1"].append(max(float(row["mass1"]), float(row["mass2"])))
        injection_values["mass2"].append(min(float(row["mass1"]), float(row["mass2"])))
        injection_values["distance"].append(float(row["distance"]))
        injection_values["psi"].append(float(row["polarization"]))
        injection_values["phase"].append(float(row["coa_phase"]))
        injection_values["geocent_time"].append(float(row["geocent_end_time"]))
        injection_values["ra"].append(float(row["longitude"]))
        injection_values["dec"].append(float(row["latitude"]))
        injection_values["spin1z"].append(float(row["spin1z"]))
        injection_values["spin2z"].append(float(row["spin2z"]))

    #injection_values = pd.DataFrame.from_dict(injection_values)
    return injection_values
    

run_names = run_dirs =  ['O3', 'O4', 'O5']

with tqdm(total=len(run_dirs)) as progress:
    for run_name, run_dir in zip(run_names, run_dirs):
        path = Path(f'{path_dir}/{run_dir}/farah/injections.xml')
        data = xml_to_dataframe(path)

        # SAve in another file format 

        from astropy.table import Table  as astro_Table

        astro_Table({"mass1"   :   data["mass1"],
                     "mass2"   :   data["mass2"],
                     "spin1z"  :   data["spin1z"],
                     "spin2z"  :   data["spin2z"],
                     "distance":   data["distance"]
                }
                   ).write(f"{outdir}/injections_xml_{run_name}.h5", overwrite=True)
        progress.update()
            


