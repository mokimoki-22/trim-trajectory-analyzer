import scanpy as sc
import pandas as pd
import numpy as np

def compute_dpt(adata, root_cell_idx):
    """
    Computes Diffusion Pseudotime (DPT) for the given AnnData object.
    Parameters:
    - adata: AnnData object
    - root_cell_idx: Index of the root cell in adata.obs_names
    """
    # Verify processing
    if 'neighbors' not in adata.uns:
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=20)
    
    # Run DPT steps
    sc.tl.diffmap(adata)
    adata.uns['iroot'] = np.where(adata.obs_names == root_cell_idx)[0][0]
    sc.tl.dpt(adata)
    
    return adata
