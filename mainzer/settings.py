import toml
import pprint


settings_list = (
# (name, type, default, comment)
('path_spectrum',                   str,    '*.mzML',           'Path to the spectrum (mzML, mzXML, csv). May contain an asterisk to match multiple files.'),
('path_base_lipids',                str,    'base_lipids.csv',  'Path to the file containing "lipid seeds".'),
('path_base_proteins',              str,    'base_proteins.csv','Path to the file containing "proteins seeds".'),
('min_highest_intensity',           float,  0.0,                'Threshold for intensities to consider [0<=x]'),
('min_mz',                          float,  0.0,                'Minimal m/z of the centroids [0<=x]'),
('max_mz',                          float,  1e20,               'Minimal m/z of the centroids [min_mz<x]'),
('min_protein_cluster_charge',      int,    1,                  'Minimal charge state for ions containing a protein.'),
('max_protein_cluster_charge',      int,    50,                 'Maximal charge state for ions containing a protein.'),
('min_lipid_mers',                  int,    1,                  'Minimal number of adducts'),
('max_lipid_mers',                  int,    4,                  'Maximal number of adducts'),
('min_free_lipid_cluster_charge',   int,    1,                  'Minimal charge state for ions containing a protein.'),
('max_free_lipid_cluster_charge',   int,    50,                 'Maximal charge state for ions containing a protein.'),
('only_heteromers',                 bool,   False,              'Consider only heteromers of proteins, without studying their monomers.'),
('min_neighbourhood_intensity',     float,  0.0,                'Minimal intensity in the neighbourhood of each ion to consider it [0<=x]'),
('min_charge_sequence_length',      int,    3,                  'Minimal number of consecutive charge states of a protein to consider it.'),
('min_total_fitted_probability',    float,  0.80,               'The minimal overall probability of isotopes matched to centroids [0<=x<=1]'),
('isotopic_coverage',               float,  0.99,               'IsoSpec probability coverage [0<=x<1] (only values close to 1, like 0.99, make sense though)'),
('isotopic_bin_size',               float,  0.1,                'IsoSpec bin size in Thomsons, [0<x]'),
('neighbourhood_thr',               float,  1.1,                'Neighbourhood buffer size in Thomsons [0<=x]'),
('max_expected_ppm_distance',       float,  20,                 'Maximal expected ppm distance between theoretical and experimental peaks between theoretical peaks and apexes of assigned centroids.'),
('underfitting_quantile',           float,  0.00,               'Single molecule underfit quantile [0<=x]'),
('min_max_intensity_threshold',     float,  100,                'Minimal maximal intensity estimate that qualifies an ion for further analysis [0<=x]'),
('chimeric_regression_fits_cnt',    int,    2,                  'Number of times to run chimeric regression.'),
('fitting_to_void_penalty',         float,  1.0,                'Penalty for fitting with theory where there is no signal [0.0<x, only used when deconvolve=true]'),
('min_chimeric_intensity_threshold',float,  100,                'In case of multiple chimeric fits, the minimal intensity of chimeric regression estimate that qualifies an ion for next iteration of the chimeric regression fitting procedure.'),
('rounding',                        int,   3,                   'Number of rounding digits. Set to -1 to disable rounding'),
('verbose',                         bool,   True,               'Verbose?'),
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

    def __repr__(self):
        return repr(self.settings)

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
