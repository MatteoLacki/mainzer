import argparse
import glob
import os
import pathlib
import sys

from mainzer.lipido import lipido_IO
from mainzer.settings import Settings

ap = argparse.ArgumentParser(description='Analyze batch input with lipido.',
                             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument('path_to_config_file',
                help="Path to configuration file.",
                type=pathlib.Path)
ap = ap.parse_args()

settings = Settings.FromTOML(ap.path_to_config_file)
os.chdir(ap.path_to_config_file.parent)


for path_spectrum in glob.glob(settings["path_spectrum"]):
    path_spectrum = pathlib.Path(path_spectrum)
    settings["path_spectrum"] = str(path_spectrum)
    output_folder = path_spectrum.stem
    lipido_IO(settings, output_folder)