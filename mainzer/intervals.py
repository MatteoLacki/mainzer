from collections import namedtuple
from ncls import FNCLS
import numpy as np
if False: # Do NOT remove this: this is to trick pyinstaller into
          # importing the necessary packages
    import scipy.spatial.transform._rotation_groups

Indices = namedtuple("IndexEdges", "query interval_db")


class IntervalQuery(object):
    """Structure for rapid building of a bi-partite graph of interval intersections.

    The intervals are open.

    The white nodes are intervals that are in the data base.
    The red nodes are intervals that are at a given moment submitted for an intersection query.
    The structure allows for efficient multiple submissions.

    Each time, the edges of the graphs are returned iff the red and white interval intersect [sparse outcome].
    This was tested agains pandas interval indices and python intervaltree and won.
    """
    def __init__(self, left, right):
        self.left  = np.array(left,  dtype=float)
        self.right = np.array(right, dtype=float)
        # ids can support indices, but this overcomplicates queries that can be simply integer based.
        self.ids = np.arange(len(left), dtype=np.int64)
        self.fncls = FNCLS(self.left, self.right, ids=self.ids)

    def interval_query(self, left, right):
        """
        Arguments:
            left (list): left ends of (open) intervals.
            right (list): right ends of (open) intervals.
        Returns:
            IndexEdges: named tuple defining edges between matching query intervals and those already in the interval database.
        """
        left = np.array(left, dtype=float)
        right = np.array(right, dtype=float)
        query_idxs, db_idxs = self.fncls.all_overlaps_both(left, right, self.ids)
        return Indices(query_idxs.astype(np.int64),
                       db_idxs.astype(np.int64))

    def interval_query_tuples(self, query_intervals):
        """
        Arguments:
            query_intervals (list): list of tuples with left and right ends of intervals.
        Returns:
            IndexEdges: named tuple defining edges between matching query intervals and those already in the interval database.
        """
        left, right = zip(*query_intervals)
        return self.interval_query(left, right)

    def point_query(self, x):
        return self.interval_query(x,x)

    def __repr__(self):
        return f"[({self.left[0]:.3f}, {self.right[0]:.3f}) .. ({self.left[-1]:.3f}, {self.right[-1]:.3f})]"


def test_IntervalQuery():
    """This can be run by typing 'pytest' from the root."""
    test_intervals = IntervalQuery([0., 2., 7.],
                                   [1., 3., 9.])
    index_edges = test_intervals.interval_query(left=[-2, 1.5],
                                                right=[-1, 10.])
    assert all(index_edges.query == [1,1]), "Wrong indices of queried intervals returned."
    assert all(index_edges.interval_db == [1,2]), "Wrong indices of intervals in the interval database returned."

    point_edges = test_intervals.point_query(x=[2.5, 6, 8])
    assert all(point_edges.query == [0,2]), "Wrong indices of queried intervals returned."
    assert all(point_edges.interval_db == [1,2]), "Wrong indices of intervals in the interval database returned." 

    query_intervals = [(-2, -1), (1.5, 10)]
    left, right = zip(*query_intervals)

    index_edges_using_tuple_query = test_intervals.interval_query_tuples(query_intervals)
    assert all(index_edges_using_tuple_query.query == index_edges.query) and \
           all(index_edges_using_tuple_query.interval_db == index_edges.interval_db), \
           "Results between 'interval_query' and 'interval_query_tuples' do not match."

#TODO: find a better place for this class.
import pandas as pd

class IonsCentroids(object):
    def __init__(self,
                 ions,
                 centroids,
                 isotopic_calculator):
        self.ions = ions
        self.centroids = centroids
        self.centroids_intervals = IntervalQuery(self.centroids.left_mz,
                                                 self.centroids.right_mz)
        self.isotopic_calculator = isotopic_calculator

    def get_isotopic_summaries(self):
        self.ions = self.isotopic_calculator.ions_summary(self.ions)

    def get_neighbourhood_intensities(self, neighbourhood_thr=1.1):
        edges = self.centroids_intervals.interval_query(self.ions.isospec_min_mz + neighbourhood_thr,
                                                        self.ions.isospec_max_mz + neighbourhood_thr)
        # different neighbourhoods might include intensity from the same centroids.
        self.ions["neighbourhood_intensity"] = self.centroids.I_sum.iloc[edges.interval_db].groupby(edges.query).sum()
        self.ions.neighbourhood_intensity.fillna(0, inplace=True)# NaN = 0 

    def assign_isotopologues_to_centroids(self, min_neighbourhood_intensity=0):
        #TODO: one might optimize this for graph purposes
        promissing_ions = self.ions[self.ions.neighbourhood_intensity > 0][["formula","charge"]]
        promissing_ions.drop_duplicates(inplace=True)
        full_envelopes = pd.concat(self.isotopic_calculator.to_frame(form, z)
                                   for form, z in promissing_ions.itertuples(index=False))
        return full_envelopes


    @staticmethod
    def from_existing_centroids(existing_ions_centroids, new_ions):
        """This is introduced to save time on unnecessary initialization."""
        return IonsCentroids(new_ions,
                             existing_ions_centroids.centroids,
                             existing_ions_centroids.isotopic_calculator)