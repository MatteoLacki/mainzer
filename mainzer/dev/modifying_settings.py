%load_ext autoreload
%autoreload 2
from mainzer.lipido import lipido_IO
from mainzer.settings import Settings

settings = Settings.FromTOML("settings.mainzer")
settings.settings["output_folder"] = "testres"

lipido_IO(settings)