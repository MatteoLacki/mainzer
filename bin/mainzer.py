# import pandas as pd
# import pathlib
# from pprint import pprint
# import json
# from datetime import datetime
# analysis_time = datetime.now().strftime('mainzer__date_%Y_%m_%d_time_%H_%M_%S')
# from nogui.forms import ask_if, ask_for

# from mainzer.read import read
# from mainzer.lipido import get_lipido_ions
# from mainzer.deconv import estimate_intensities


# print()

# logo = """                                                                                                    
#               &@@                                                                                   
#              @@ ,@#                                                         .@,                     
#            *@@*   @%                                                        &@&                     
#           .@@@    /@                                 @@@@@@@@@@@@           @@@                     
#           @@@.    .@*      *%#,                    /@         *@%          *@@@                     
#          %@@@     ,@.   ,@@,   &@*                    .%&   %@%            %@@@                     
#          @@@,     &@    @/     @@                  &@@  (@@&               @@@@                     
#         &@@@      @#         @@              ,&@@@                         @@@&                     
#         @@@/     @@   .@@@@&            &@@@/     @@     @@                @@@%          .&@@@@%    
#        %@@@    (@(        /@       @@@@(           @@    @&           *@@% @@@#      @@@%      @@   
#        @@@&   ,.          @@     (@@@              @@   .@%        #@@     @@@*    @@/         @&   
#       ,@@@,               @@    @*@@@             &@(   ,@#      &@/       @@@.  ,@/*%       (@(    
#       @@@@               .@#     .@@@            @@,    *@(     %@       #&@@@   @@        %@&      
#       @@@@               %@       @@@         .@@&      *@(      @@@%%@@   @@@   @@     &@@         
#      (@@@*               &@      /@@@     ,&@@@         *@/                @@@      *(*             
#      @@@@                        &@@@   @#              *@,                @@@                      
#      @@@@                        @@@&                                       @@                      
#     .@@@&          @@@@@@@@@&    @@@@                                                               
#     %@@@*      .@@@@@@#.   &@@(  @@@@                                                               
#     @@@@     #@@@@@          &@# @@@@                                                               
#     @@@@   &@@@@(             /  @@@@                  
#     @@@@ *@@@@*                  @@@&                 
#    /@@@*@@@@%                    @@@&                 
#    &@@@@@@@                      @@@%                
#    @@@@@@(                       @@@(                
#    @@@@@                         @@@#                
#   ,@@@&                           @@@,               
#   @@@(                                               
#   @@@                                                                                               
#                                                                                                    """
# print(logo)
# print("This program means business!")
# print()

# # Make it accept CLI arguments!

# path_ions = pathlib.Path(input("Paste in the path to 'ions.csv': "))
# # path_molecules = pathlib.Path("/home/matteo/Projects/och_Kallol/unlipid/data/test/molecules.csv")
# assert path_ions.exists(), "The csv with ions is not available."
# ions = pd.read_csv(path_ions)
# assert all(colname in ions.columns for colname in ('name','charge','formula'))

# path_spectrum = pathlib.Path(input("Paste in the path to the spectrum: "))
# # path_spectrum = pathlib.Path("/home/matteo/Projects/och_Kallol/unlipid/data/07232020_Resolution50000_64bit.mzML")
# assert path_spectrum.exists(), "This spectrum file is not available."
# mz, intensity = read(path_spectrum)


# settings = {}
# settings["path_molecules"] = str(path_molecules)
# settings["path_spectrum"] = str(path_spectrum)
# output_folder = pathlib.Path(input("Ouput folder: ")).expanduser()
# output_folder.mkdir(parents=True, exist_ok=True)
# (output_folder/analysis_time).mkdir(parents=True, exist_ok=True)
# settings["output_folder"] = str(output_folder)

# settings["isotopic_coverage"] = ask_for("IsoSpec probability coverage [0<x<1]:", .95, float)
# settings["isotopic_bin_size"] = ask_for("IsoSpec bin size in Thomsons [0<x]:", .1, float)
# settings["neighbourhood_thr"] = ask_for("Neighbourhood buffer size in Thomsons [0<x]:", 1.1, float)
# settings["underfitting_quantile"] = ask_for("Single molecule underfit quantile [0<x]:", .05, float)
# settings["deconvolve"] = ask_if("Deconvolve? [Y/n]")
# if settings["deconvolve"]:
#     settings["fitting_to_void_penalty"] = ask_for("Penalty for fitting with theory where there is no signal [0<x]:", 1.0, float)
# settings["verbose"] = ask_if("Verbose? [Y/n]")

# print()
# print("Running Mainzer with:")
# pprint(settings)
# print()
# print("It's business time!")
# print("Estimating intenisities")
# ions, timings = estimate_intensities(mz, intensity, ions, verbose_output=False, **settings)

# column_order = ["name"]
# if "deconvolved_intensity" in ions.columns:
#     ions = ions.sort_values(['charge','deconvolved_intensity'], ascending=[True, False])
#     column_order.append("deconvolved_intensity")
# else:
#     ions = ions.sort_values(['charge','maximal_intensity'], ascending=[True, False])
    
# column_order.extend(["maximal_intensity",
#                      "proximity_intensity",
#                      "neighbourhood_intensity",
#                      "isospec_final_coverage",
#                      "isospec_prob_with_signal",
#                      "isospec_prob_without_signal",
#                      "isospec_peaks_count",
#                      "min_isospec_mz",
#                      "max_isospec_mz"])

# ions = ions[column_order]
# final_folder = output_folder/analysis_time

# print("Saving results")
# ions.to_csv(final_folder/"ions.csv")

# with open(final_folder/"timings.json", "w") as jsonfile:
#     json.dump(timings, jsonfile, indent=4)
# with open(final_folder/"settings.json", "w") as jsonfile:
#     json.dump(settings, jsonfile, indent=4)

# print("Thank you for letting Mainzer do its job!")
