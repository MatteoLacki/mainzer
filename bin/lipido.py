import pandas as pd
import pathlib
from pprint import pprint
import json
from datetime import datetime
analysis_time = datetime.now().strftime('lipido__date_%Y_%m_%d_time_%H_%M_%S')
from nogui.forms import ask_if, ask_for

import IsoSpecPy # Not used, here, but to make sure pyinstaller pulls it in

from mainzer.read import read
from mainzer.lipido import get_lipido_ions
from mainzer.deconv import estimate_intensities
from mainzer.settings import Settings


print()

logo = """                                                                                                    
              &@@                                                                                   
             @@ ,@#                                                         .@,                     
           *@@*   @%                                                        &@&                     
          .@@@    /@                                 @@@@@@@@@@@@           @@@                     
          @@@.    .@*      *%#,                    /@         *@%          *@@@                     
         %@@@     ,@.   ,@@,   &@*                    .%&   %@%            %@@@                     
         @@@,     &@    @/     @@                  &@@  (@@&               @@@@                     
        &@@@      @#         @@              ,&@@@                         @@@&                     
        @@@/     @@   .@@@@&            &@@@/     @@     @@                @@@%          .&@@@@%    
       %@@@    (@(        /@       @@@@(           @@    @&           *@@% @@@#      @@@%      @@   
       @@@&   ,.          @@     (@@@              @@   .@%        #@@     @@@*    @@/         @&   
      ,@@@,               @@    @*@@@             &@(   ,@#      &@/       @@@.  ,@/*%       (@(    
      @@@@               .@#     .@@@            @@,    *@(     %@       #&@@@   @@        %@&      
      @@@@               %@       @@@         .@@&      *@(      @@@%%@@   @@@   @@     &@@         
     (@@@*               &@      /@@@     ,&@@@         *@/                @@@      *(*             
     @@@@                        &@@@   @#              *@,                @@@                      
     @@@@                        @@@&                                       @@                      
    .@@@&          @@@@@@@@@&    @@@@                                                               
    %@@@*      .@@@@@@#.   &@@(  @@@@                                                               
    @@@@     #@@@@@          &@# @@@@                                                               
    @@@@   &@@@@(             /  @@@@                  
    @@@@ *@@@@*                  @@@&                 
   /@@@*@@@@%                    @@@&                 
   &@@@@@@@                      @@@%                
   @@@@@@(                       @@@(                
   @@@@@                         @@@#                
  ,@@@&                           @@@,               
  @@@(                                               
  @@@                                                                                               
                                                                                                   """
print(logo)
print("This program means business!")
print()

# Make it accept CLI arguments!
'''
path_molecules = pathlib.Path(input("Paste in the path to 'molecules.csv': "))
# path_molecules = pathlib.Path("/home/matteo/Projects/och_Kallol/unlipid/data/test/molecules.csv")
assert path_molecules.exists(), "The csv with molecules is not available."

path_spectrum = pathlib.Path(input("Paste in the path to the spectrum: "))
# path_spectrum = pathlib.Path("/home/matteo/Projects/och_Kallol/unlipid/data/07232020_Resolution50000_64bit.mzML")
assert path_spectrum.exists(), "This spectrum file is not available."
'''

settings = Settings.FromConsole()
'''
settings["path_molecules"] = str(path_molecules)
settings["path_spectrum"] = str(path_spectrum)
settings["output_folder"] = str(output_folder)
settings["min_intensity"] = ask_for("Threshold for intensities to consider [0<x]:", 0.0, float)
settings["max_lipid_mers"] = ask_for("Maximal number of adducts:", 4)
settings["min_protein_charge"] = ask_for("Minimal protein charge state:", 1)
settings["max_protein_charge"] = ask_for("Maximal protein charge state:", 10)
settings["min_lipid_charge"] = ask_for("Minimal charge on a free lipid cluster:", 1)
settings["max_lipid_charge"] = ask_for("Maximal charge on a free lipid cluster:", 10)
settings["isotopic_coverage"] = ask_for("IsoSpec probability coverage [0<x<1]:", .95, float)
settings["isotopic_bin_size"] = ask_for("IsoSpec bin size in Thomsons [0<x]:", .1, float)
settings["neighbourhood_thr"] = ask_for("Neighbourhood buffer size in Thomsons [0<x]:", 1.1, float)
settings["underfitting_quantile"] = ask_for("Single molecule underfit quantile [0<x]:", .05, float)
settings["deconvolve"] = ask_if("Deconvolve? [Y/n]")
if settings["deconvolve"]:
    settings["fitting_to_void_penalty"] = ask_for("Penalty for fitting with theory where there is no signal [0<x]:", 1.0, float)
settings["verbose"] = ask_if("Verbose? [Y/n]")
'''
molecules = pd.read_csv(settings['path_molecules'])
mz, intensity = read(settings['path_spectrum'])
output_folder = pathlib.Path(input("Ouput folder: ")).expanduser()
output_folder.mkdir(parents=True, exist_ok=True)
(output_folder/analysis_time).mkdir(parents=True, exist_ok=True)

mz, intensity = mz[intensity > settings["min_intensity"]], intensity[intensity > settings["min_intensity"]]

print()
print("Running Lipido with:")
pprint(settings)
print()
print("It's business time!")

print("Getting ions")
ions = get_lipido_ions(molecules, **settings)

print("Estimating intenisities")
ions, centroids, timings  = estimate_intensities(mz, intensity, ions, verbose_output=False, **settings)

column_order = ["name"]
if "deconvolved_intensity" in ions.columns:
    ions = ions.sort_values(['charge','deconvolved_intensity'], ascending=[True, False])
    column_order.append("deconvolved_intensity")
else:
    ions = ions.sort_values(['charge','maximal_intensity'], ascending=[True, False])
    
column_order.extend(["maximal_intensity",
                     "proximity_intensity",
                     "neighbourhood_intensity",
                     "isospec_final_coverage",
                     "isospec_prob_with_signal",
                     "isospec_prob_without_signal",
                     "touched_centroids",
                     "isospec_peaks_count",
                     "min_isospec_mz",
                     "max_isospec_mz"])

ions = ions[column_order]
protein_names = "|".join(molecules.query("group == 'protein'").name)
(_, proteins), (_, lipid_clusters) = ions.groupby(ions.name.str.contains(protein_names), group_keys=False)

final_folder = output_folder/analysis_time
# ions.to_csv(/"ions.csv")

print("Saving results")
proteins.to_csv(final_folder/"proteins.csv")
lipid_clusters.to_csv(final_folder/"lipid_clusters.csv")
centroids.df.to_csv(final_folder/"centroids.csv")

with open(final_folder/"timings.json", "w") as jsonfile:
    json.dump(timings, jsonfile, indent=4)
with open(final_folder/"settings.json", "w") as jsonfile:
    json.dump(settings, jsonfile, indent=4)

print("Thank you for letting Lipido do its job!")
