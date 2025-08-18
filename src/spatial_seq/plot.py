import pandas as pd
import numpy as np
from glasbey import create_palette
from typing import Iterable
from scanpy import AnnData


def color_gen(groups: pd.Series | Iterable | np.array, custom_index=None):
    cs = create_palette(palette_size=len(groups.unique()))
    if custom_index is not None:
        return pd.Series(cs, index=custom_index)
    else:
        return pd.Series(cs, index=groups.unique())


def order_obs(adata: AnnData, col: str, order: Iterable[str]):
    adata.obs[col] = pd.Categorical(adata.obs[col], categories=order, ordered=True)
    return
