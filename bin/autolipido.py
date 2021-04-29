import sys
import os
import pandas as pd
from pathlib import Path
from pprint import pprint
import json
import toml
import glob
from datetime import datetime
analysis_time = datetime.now().strftime('lipido__date_%Y_%m_%d_time_%H_%M_%S')
from nogui.forms import ask_if, ask_for

from mainzer.read import read
from mainzer.lipido import lipido_main
from mainzer.deconv import estimate_intensities
from mainzer.settings import Settings


def usage():
    print('''
Usage:
    python autolipido.py <path_to_config_file>
''')

if len(sys.argv) < 2:
    usage()
    sys.exit(1)


config_path = Path(sys.argv[1])

settings = Settings.FromTOML(config_path)

os.chdir(config_path.parent)


for path_spectrum in glob.glob(settings["path_spectrum"]):
    path_spectrum = Path(path_spectrum)
    settings["path_spectrum"] = str(path_spectrum)
    settings['output_folder'] = str(Path(path_spectrum.stem) / "output")
    lipido_main(settings)
