import pandas as pd
import pathlib
from pprint import pprint
import json

import IsoSpecPy # Not used, here, but to make sure pyinstaller pulls it in

from mainzer.read import read
from mainzer.lipido import lipido_main
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

settings = Settings.FromConsole()
settings['output_folder'] = str(pathlib.Path(input("Ouput folder: ")).expanduser())
lipido_main(settings)
