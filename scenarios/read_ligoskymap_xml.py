import pandas as pd
from gwpy.table import Table as gwpy_Table
from astropy.table import Table  as astro_Table

# Read Ligo xml file generate by ligo.skymap


    
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

fara_xml_dir = "/home/weizmann.kiendrebeogo/OBSERVING_SCENARIOS/observing-scenarios-2022/XML/Farah"

leo_xml_dir =  "/home/weizmann.kiendrebeogo/OBSERVING_SCENARIOS/observing-scenarios-2022/XML/Leo"

file_dir =[fara_xml_dir, leo_xml_dir]   

runs = ['O3', 'O4', 'O5']
pops = ['BNS','NSBH','BBH']

for file in file_dir:
    for run in runs:
        for pop in pops:
            if file == fara_xml_dir:
                pop_dir = f"{file}/{run}/farah_{pop.lower()}" 
            else:
                pop_dir = f"{file}/{run}/{pop.lower()}_astro"
                
            xml_path = f"{pop_dir}/injections.xml"
            data = xml_to_dataframe(xml_path)
            
            astro_Table({
                "mass1"   :   data["mass1"],
                "mass2"   :   data["mass2"],
                "spin1z"  :   data["spin1z"],
                "spin2z"  :   data["spin2z"],
                "distance":   data["distance"]
                }
           ).write(f"{pop_dir}/injections.h5", overwrite=True)




# SAve in another file format 




 
