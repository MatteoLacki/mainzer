def settings2toml(path_spectrum = '*.mzML',
                  path_molecules = 'molecules.csv',
                  max_lipid_mers = 4,
                  min_protein_charge = 1,
                  max_protein_charge = 10,
                  min_lipid_charge = 1,
                  max_lipid_charge = 10,
                  isotopic_coverage = 0.99,
                  isotopic_bin_size = 0.1,
                  neighbourhood_thr = 1.1,
                  underfitting_quantile = 0.05,
                  deconvolve = True,
                  fitting_to_void_penalty = 1.0,
                  verbose = True,
                  **kwds):
    """Save typical settings to a TOML template for 'config.mainzer'."""
    deconvolve = 'true' if deconvolve else 'false'
    verbose = 'true' if verbose else 'false'
    toml_str = f"""
    # Every line starting with a # is a comment

    # Star allowed only here, will process all the matching files
    # Path should be absolute, or relative to the .mainzer file
    # Path to the spectrum (mzML, mzXML, csv)
    path_spectrum = {path_spectrum}

    # Path to the file containing "ion seeds": molecules that give rise to all the ions we search for
    path_molecules = {path_molecules}

    # Maximal number of adducts
    max_lipid_mers = {max_lipid_mers}

    # Minimal protein charge state
    min_protein_charge = {min_protein_charge}

    # Maximal protein charge state
    max_protein_charge = {max_protein_charge}

    # Minimal charge on a free lipid cluster
    min_lipid_charge = {min_lipid_charge}

    # Maximal charge on a free lipid cluster
    max_lipid_charge = {max_lipid_charge}

    # IsoSpec probability coverage [0<x<1] (only values close to 1, like 0.99, make sense though)
    isotopic_coverage = {isotopic_coverage}

    # IsoSpec bin size in Thomsons [0<x]
    isotopic_bin_size = {isotopic_bin_size}

    # Neighbourhood buffer size in Thomsons [0<x]
    neighbourhood_thr = {neighbourhood_thr}

    # Single molecule underfit quantile [0<x]
    underfitting_quantile = {underfitting_quantile}

    # Deconvolve? [true/false]
    deconvolve = {deconvolve}

    # Penalty for fitting with theory where there is no signal [0<x, only used when deconvolve=true]
    fitting_to_void_penalty = {fitting_to_void_penalty}

    # Verbose? [true/false]
    verbose = {verbose}
    """
    return toml_str
