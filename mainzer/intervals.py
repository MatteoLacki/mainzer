from ncls import FNCLS
import numpy as np


class IntervalQuery(object):
    """Structure for rapid building of a bi-partite graph of interval intersections.

    The intervals are open.

    The white nodes are intervals that are in the data base.
    The red nodes are intervals that are at a given moment submitted for an intersection query.
    The structure allows for efficient multiple submissions.

    Each time, the edges of the graphs are returned iff the red and white interval intersect [sparse outcome].
    This was tested agains pandas interval indices and python intervaltree and won.
    """
    def __init__(self, left, right, index=None):
        self.left  = np.array(left,  dtype=float)
        self.right = np.array(right, dtype=float)
        index = np.arange(len(left)) if index is None else np.array(index)
        self.fncls = FNCLS(self.left, self.right, index)

    def interval_query(self, left, right, index=None):
        left = np.array(left, dtype=float)
        right = np.array(right, dtype=float)
        if index == None:
            index = np.arange(len(left), dtype=np.int64)
        else:
            index = index.astype(np.int64)
        query_idxs, db_idxs = self.fncls.all_overlaps_both(left, right, index)
        return query_idxs.astype(np.int64), db_idxs.astype(np.int64)

    def point_query(self, x, index=None):
        return self.interval_query(x,x,index)

    def __repr__(self):
        return f"[({self.left[0]:.3f}, {self.right[0]:.3f}) .. ({self.left[-1]:.3f}, {self.right[-1]:.3f})]"



# if __name__ == "__main__":
#     test_intervals = IntervalQuery([0., 2., 7.],
#                                    [1., 3., 9.])
#     print(test_intervals.interval_query([1.5], [10.])
