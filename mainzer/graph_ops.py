import collections
import networkx as nx
import numpy as np
import pandas as pd


def bipartiteGraph2regressionProblem(connected_component,
                                     merge_zeros=True,
                                     normalize_X=True):
    Y = pd.Series(connected_component.intensities(), name="intensity")
    prob_out = pd.Series(connected_component.unmatched_probability())
    X_np = nx.algorithms.bipartite.matrix.biadjacency_matrix(G=connected_component,
                                                             weight="prob",
                                                             row_order=Y.index,
                                                             column_order=prob_out.index).toarray()
    X = pd.DataFrame(X_np, columns=prob_out.index, index=Y.index)
    if merge_zeros:
        X_zero = pd.DataFrame([prob_out.values.T], columns=prob_out.index, index=[-1] )
        Y_zero = pd.Series([0], index=[-1], name="intensity")
    else:
        X_zero = pd.DataFrame( np.diag(prob_out), columns=prob_out.index, index=[-1]*len(prob_out) )
        Y_zero = pd.Series(np.zeros(shape=len(X_zero)), index=[-1]*len(X_zero), name="intensity")
    Y_full = pd.concat([Y, Y_zero])
    X_full = pd.concat([X, X_zero])
    if normalize_X:
        X_full = X_full.dot( np.diag( 1/X_full.sum(axis=0) ) )
        X_full.columns = prob_out.index
    return X_full, Y_full



class RegressionGraph(nx.Graph):
    def connected_components(self):
        for cc in nx.connected_components(self):
            yield self.subgraph(cc)

    def count_connected_components(self):
        i = 0
        for cc in nx.connected_components(self):
            i += 1
        return i

    def count_nodes_in_subproblems(self):
        problem_size = lambda x: "interval_cnt" if isinstance(x, int) else "charged_formula_cnt"
        return pd.DataFrame(collections.Counter(problem_size(e) for e in cc) for cc in nx.connected_components(self))

    def get_problem_sizes(self):
        return pd.DataFrame((X.shape for X,Y in self.iter_regression_problems()),
                            columns=("interval_cnt","charged_formula_cnt"))

    def iter_convoluted_ions(self):
        for cc in nx.connected_components(self):
            yield {node for node in cc if isinstance(node, tuple)} 

    def plot_problem_sizes(self, show=True):
        import matplotlib.pyplot as plt
        complexity = self.get_problem_sizes() 
        plt.scatter(complexity.interval_cnt, complexity.charged_formula_cnt)
        plt.xlabel("Number of intevals")
        plt.ylabel("Number of explanations")
        plt.xscale('log')
        plt.yscale('log')
        if show:
            plt.show()

    def iter_regression_problems(self, merge_zeros=True, normalize_X=True):
        for cc in self.connected_components():
            yield bipartiteGraph2regressionProblem(cc, merge_zeros=merge_zeros, normalize_X=normalize_X)

    def intensities(self):
        return nx.get_node_attributes(self, "intensity")

    def unmatched_probability(self):
        return nx.get_node_attributes(self, "prob_out")

    def iter_chimeric_groups(self):
        for idx, ci in enumerate(self.iter_convoluted_ions()):
            formulas, charges = zip(*ci)    
            yield pd.DataFrame({"formula": formulas,
                                "charge": charges,
                                "chimeric_group": idx})