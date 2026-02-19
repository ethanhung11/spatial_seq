# uv run python src/scripts/create_cloupe.py &> ./outs/cloupe.out

from pathlib import Path
import scanpy as sc
from utils import create_cloupe
adata = sc.read_h5ad(Path("../../data") / "processed" / "single_cell" / "combined" / "eWAT_Male.h5ad")

create_cloupe(
    adata,
    layer = "counts",
    obs = ["Identifier", "Dataset", "Diet", "Age", "Sex", "Depot", "celltype"],
    obsm = [
        "UMAP_INT_none",
        "UMAP_INT_scvi_hvg-Identifier",
        "LocalMAP_INT_scvi_hvg-Identifier",
    ],
    output = Path("data") / "processed" / "single_cell" / "combined" / "eWAT_Male_adipose_atlas",
    new_barcodes = None,
)