import toml
import pprint


settings_list = (
# (name, type, default, comment)
('path_spectrum',           str,    '*.mzML',           'Path to the spectrum (mzML, mzXML, csv). May contain an asterisk to match multiple files.'),
('path_molecules',          str,    'molecules.csv',    'Path to the file containing "ion seeds": molecules that give rise to all the ions we search for'),
('min_intensity',           float,  0.0,                'Threshold for intensities to consider [0.0<=x]'),
('max_protein_mers',        int,    2,                  'Maximal oligomeric state of protein'),
('max_lipid_mers',          int,    4,                  'Maximal number of adducts'),
('min_protein_charge',      int,    1,                  'Minimal protein charge state'),
('max_protein_charge',      int,    10,                 'Maximal protein charge state'),
('min_lipid_charge',        int,    1,                  'Minimal charge on a free lipid cluster'),
('max_lipid_charge',        int,    10,                 'Maximal charge on a free lipid cluster'),
('isotopic_coverage',       float,  0.99,               'IsoSpec probability coverage [0.0<x<1.0] (only values close to 1, like 0.99, make sense though)'),
('isotopic_bin_size',       float,  0.1,                'IsoSpec bin size in Thomsons, [0.0<x]'),
('neighbourhood_thr',       float,  1.1,                'Neighbourhood buffer size in Thomsons [0.0<x]'),
('underfitting_quantile',   float,  0.05,               'Single molecule underfit quantile [0.0<x]'),
('deconvolve',              bool,   True,               'Deconvolve?'),
('fitting_to_void_penalty', float,  1.0,                'Penalty for fitting with theory where there is no signal [0.0<x, only used when deconvolve=true]'),
('verbose',                 bool,   True,               'Verbose?'),
)

setting_types = { name : t for name, t, _, _ in settings_list }
setting_descr = { name : desc for name, _, _, desc in settings_list }

class Settings:
    def __init__(self):
        self.settings = { name: value for name, _, value, _ in settings_list }

    @staticmethod
    def FromTOML(path_or_filelike):
        try:
            with open(path_or_filelike) as f:
                toml_str = f.read()
        except TypeError:
            toml_str = path_or_filelike.read()

        result = Settings()
        result.settings.update(toml.loads(toml_str))
        return result

    @staticmethod
    def FromConsole():
        result = Settings()
        def bool_conv(s):
            if s.lower() in "y t true yes".split():
                return True
            if s.lower() in "n f false no".split():
                return False
            raise ValueError("Response " + s + " not understood.")

        for name, type_, default, descr in settings_list:
            conv = type_
            if type_ == bool:
                descr += " [Y/n] >" if default else " [N/y] > "
                conv = bool_conv
            else:
                descr = f"{descr} (default: {default}) > "

            while True:
                try:
                    print(descr, end='')
                    i = input()
                    val = default if i == '' else conv(i)
                    result[name] = val
                    break
                except ValueError as e:
                    print(str(e))

        return result


    def __getitem__(self, key):
        return self.settings[key]

    def __setitem__(self, key, val):
        self.settings[key] = val

    def __str__(self):
        return pprint.pformat(self.settings)

    def print_summary(self):
        print(str(self))

    def toml_str(self):
        result = ["# Every line starting with a # is a comment\n"]

        for name, t, _, descr in settings_list:
            if t == bool:
                descr += " [true/false]"
                val = "true" if self[name] else "false"
            elif t in (int, float):
                val = str(self[name])
            elif t == str:
                val = "'" + self[name] + "'"
            else:
                raise RuntimeError("Unknown type")

            result.append("\n# ")
            result.append(descr)
            result.append("\n")
            result.append(name)
            result.append(" = ")
            result.append(val)
            result.append("\n")

        return ''.join(result)

    def save_toml(self, path_or_filelike):
        toml_s = self.toml_str()
        try:
            with open(path_or_filelike, "w") as f:
                f.write(toml_s)
        except TypeError:
            path_or_filelike.write(toml_s)
            
if __name__ == "__main__":
    print(Settings.FromConsole().toml_str())
