%load_ext autoreload
%autoreload 2
import mainzer.isotope_ops

isotopic_calculator = mainzer.isotope_ops.IsotopicCalculator(
    coverage=.999,
    bin_size=.2,
)

isotopic_calculator("C100H202").plot()
isotopic_calculator.masses_probs("C100H202")
