"""This is place for some more general fragments."""
import pandas as pd
from typing import Iterable


def get_roepstorrf_scheme_fragments(
    proteins: pd.DataFrame,
    fragments: str="by",
) -> Iterable:
    """Generate fragments.

    Arguments:
        proteins (pd.DataFrame): A data frame with columns "name" and "sequence".
        fragments (str): Types of fragments to generate.
    """

    for fragment_type in fragments:
        pass
    pass
