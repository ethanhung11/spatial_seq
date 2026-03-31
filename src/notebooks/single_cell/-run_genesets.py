# pixi run -e main python -u ./src/notebooks/single_cell/-run_genesets.py &> outs/single_cell/run_genesets.out

# base
import sys
from pathlib import Path
import session_info
import scanpy as sc

# custom
sys.path.insert(0, 'src')
# from single_cell.R import *
# from single_cell.preprocess import *
from single_cell.plot import *
from single_cell.analysis import *
from utils import *

CORES = 10
DATADIR = Path("data") / "processed" / "single_cell" / "combined"
REFDIR = Path("references")
SUBSETS_DIR = DATADIR / "subsets"
CELLTYPE_ADDON = "fibro"

if __name__ == "__main__":
    session_info.show()
    START = stopwatch(mode=2)
    adata_fibro = sc.read_h5ad(SUBSETS_DIR / (CELLTYPE_ADDON + '.h5ad'))

    # import genesets
    import gseapy as gp
    mouse_hallmark = gp.Msigdb().get_gmt(category="mh.all", dbver="2026.1.Mm")
    genesets = {
        "Mouse_HM" : mouse_hallmark
    }

    for gs_name in genesets:
        for method in tqdm(["Seurat", "UCell", "Decoupler"]):
            if method == "Decoupler":
                for submethod in tqdm(["ULM", "AUCell", "FGSEA"]):
                    module_score(adata_fibro, genesets[gs_name], gs_name, method, submethod)
            else:
                module_score(adata_fibro, genesets[gs_name], gs_name, method)
    
    adata_fibro.write(SUBSETS_DIR / (CELLTYPE_ADDON + '-GS.h5ad'))

    stopwatch("SCRIPT", START)