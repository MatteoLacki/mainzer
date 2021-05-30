%load_ext autoreload
%autoreload 2

import numpy as np
import json
from mainzer.read import read_spectrum
from mainzer.signal_ops import cluster_spectrum
from mainzer.intervals import IntervalQuery
from mainzer.centroiding import QueryCentroids

# fixing clustering
with open("mainzer/dev/settings.json") as f:
    settings = json.load(f)
mz, intensity = read_spectrum(settings['path_spectrum'])
min_intensity_threshold = 100

clusters_df = cluster_spectrum(mz, intensity)
# this is simple: filerting on `max_intensity` in `clusters_df`
filtered_clusters_df = clusters_df[clusters_df.I_max >= min_intensity_threshold].copy()

IQ = IntervalQuery(filtered_clusters_df.left_mz,
                   filtered_clusters_df.right_mz)

IQ.point_query([100, 400, 1000, 398.73])
IQ.interval_query([100, 400, 1000, 398.73], [101, 401, 1001, 398.83])

intervals = np.array([(100,101), (400,401), (1000,1001), (398.73, 398.83)])
res = IQ.interval_query_tuples(intervals)
intervals[res.query]
filtered_clusters_df.iloc[res.interval_db]


res = isotopic_envelopes.envelopes_summary()
ions_df.merge(res)
isotopic_envelopes.ions_summary()


res["charge"] = ions_df.charge

ions_centroids = IonsCentroids(protein_ions,
                               filtered_clusters_df,
                               isotopic_calculator)
ions_centroids.get_isotopic_summaries()
ions_centroids.get_neighbourhood_intensities(neighbourhood_thr)
ions_centroids.assign_isotopologues_to_centroids()

