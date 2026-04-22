# pixi run -e main python -u ./src/notebooks/single_cell/1.2-fibroblasts_int_2.py &> outs/single_cell/fibro_integration.out

# base
import sys
from pathlib import Path
import session_info
import scanpy as sc

# custom
sys.path.insert(0, 'src')
from single_cell.preprocess import *
from utils import *

CORES = 10
DATADIR = Path("data") / "processed" / "single_cell" / "combined"
REFDIR = Path("references")
SUBSETS_DIR = DATADIR / "subsets"
CELLTYPE_ADDON = "fibro"

if __name__ == "__main__":

    session_info.show()
    START = stopwatch(mode=-1)

    task = f"LOAD DATA"
    TASK_START = stopwatch(task, START, mode=0)

    adata_fibro = sc.read_h5ad(SUBSETS_DIR / "fibro.h5ad")
    Normalize(adata_fibro)
    adata_fibro.X = adata_fibro.layers["normalized"].copy()

    stopwatch(task, TASK_START, mode=1)

    for integration_method in ["harmony"]:

        for batch_column in ["Dataset", "Identifier"]:

            pca_key = f"global_PCA-hvg"
            int_key = f"INT_{integration_method}-{pca_key}-{batch_column}"

            task = f"INTEGRATION `{integration_method}` with `{pca_key}` genes over label `{batch_column}`"
            TASK_START = stopwatch(task, START, mode=0)

            # Sort by batch to be contiguous (scanorama is picky or whateva)
            idx = adata_fibro.obs.sort_values(batch_column).index
            adata_fibro = adata_fibro[idx]

            adata_fibro = Integrate(
                adata_fibro,
                batch_column,
                kind=integration_method,
                pca_key=pca_key,
                integration_key=int_key,
            )

            print(adata_fibro.obsm)
            assert int_key in list(adata_fibro.obsm)

            Visualize(adata_fibro, key=int_key, obsm=int_key, localmap=False, show=False)
            stopwatch(task, TASK_START, mode=1)

        adata_fibro.write(SUBSETS_DIR / f"{CELLTYPE_ADDON}_reintegrated_2.h5ad")

        for N in [500,2000]:
            
            pca_key = f"{CELLTYPE_ADDON}_PCA-hvg-{N}"
            FindVariableGenes(adata_fibro, "seurat_v3", n_features=500)
            PCA(adata_fibro, key=pca_key)

            for batch_column in ["Dataset", "Identifier"]:
                task = f"INTEGRATION Method `{integration_method}` with `{pca_key}` genes over label `{batch_column}`"
                TASK_START = stopwatch(task, START, mode=0)

                # Sort by batch to be contiguous (scanorama is picky or whateva)
                idx = adata_fibro.obs.sort_values(batch_column).index
                adata_fibro = adata_fibro[idx]

                # Run integration
                int_key = f"INT_{integration_method}-{pca_key}-{batch_column}"
                Integrate(
                    adata_fibro,
                    batch_column,
                    pca_key=pca_key,
                    kind=integration_method,
                    integration_key=int_key,
                    verbose=True,
                )
                Visualize(adata_fibro, key=int_key, obsm=int_key, localmap=False, show=False)
                stopwatch(task, TASK_START, mode=1)

        adata_fibro.write(SUBSETS_DIR / f"{CELLTYPE_ADDON}_reintegrated_2.h5ad")

print("\n\nScript Complete!")