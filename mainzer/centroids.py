import numpy as np

from .signal_ops import cluster_spectrum
from .intervals import IntervalQuery


class Centroids(object):
    def __init__(self, mz, intensity):
        self.mz = np.array(mz)
        self.intensity = np.array(intensity)
        self.df = cluster_spectrum(mz, intensity).drop(['left_idx','right_idx'], axis=1)
        self.df.index.name = 'cluster'
        # it is important to add the self.df.index:
        # otherwise the df.loc won't work, if we somehow filtered df. 
        self.index = IntervalQuery(self.df.left_mz, self.df.right_mz, self.df.index)

    def __repr__(self):
        return "Clusters:\n" + repr(self.df)

    def interval_query(self, left_mz, right_mz, use_query_index=True, add_cluster_index=True, add_query_mzs=True):
        query_idx, df_idx = self.index.interval_query(left_mz, right_mz)
        res = self.df.loc[df_idx].reset_index()
        if use_query_index:
            try:
                assert np.all(left_mz.index == right_mz.index), "Indices are different where they should be the same!!!"
                res.index = left_mz.iloc[query_idx].index
            except AttributeError:
                res.index = query_idx
                add_query_mzs = False
        else:
            res.index = query_idx
        if add_query_mzs:
            res["query_left_mz"] = left_mz.iloc[query_idx].values
            res["query_right_mz"] = right_mz.iloc[query_idx].values
        if add_cluster_index:
            res = res.set_index('cluster', append=True)
        return res

    def point_query(self, mz, use_query_index=True, use_df_index=False, add_query_mzs=True):
        res = self.interval_query(mz, mz, False, False, False)
        if add_query_mzs:
            found_ms = mz.iloc[res.index]
            res.index = found_ms.index 
            res["query_mz"] = found_ms
        if use_df_index:
            res = res.set_index('cluster', append=True)
        return res 
