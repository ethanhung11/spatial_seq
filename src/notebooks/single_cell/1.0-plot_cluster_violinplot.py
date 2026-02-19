# uv run python src/notebooks/single_cell/1.0-plot_cluster_violinplot.py &> ./outs/cluster_violin_plots.out

# base
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad
from single_cell.plot import plot_cluster_violinplot
from utils import stopwatch
import os

import logging

logging.basicConfig(level=logging.ERROR)

CORES = 10
DATADIR = Path("data")
PLOTDIR = Path("plots")
REFDIR = Path("references")
MAIN_DIR = DATADIR / "processed" / "single_cell" / "combined"
SUBSETS_DIR = MAIN_DIR / "subsets"

METADATA = ["Diet", "Age", "Depot", "Sex"]
DOUBLETMETHODS = ["scDblFinder", "DoubletFinder", "doubletdetection", "scrublet"]
INT_KEY = "INT_harmony-Identifier"

if __name__ == "__main__":
    # load
    adata_macro = ad.read_zarr(SUBSETS_DIR / "macro")
    adata_fibro = ad.read_zarr(SUBSETS_DIR / "fibro")

    mac_paper_markers = {
        "Arg1 FALC Macs": [
            "Arg1",
            "Thbs1",
            "Lyz1",
            "Cd5l",
            "Mertk",
            "Timd4",
            "Fn1",
            "Clec4d",
            "Marco",
            "Itpr1",
            "Slpi",
            "Alox15",
            "Cd38",
        ],
        "Lipo VAMs": ["Maf", "Lpl", "Mrc1", "H2-Ab1", "Lipa", "Retnla", "Ly6e"],
        "Septal Macs": ["Mrc1", "Maf", "Lyve1", "Cd209f", "Cd74", "Mgl2"],
        "Classical Mono": ["Ccr2", "Ly6c1", "Ly6c2", "Cd226", "Egfr"],
        "Non-classical Mono": ["Plac8"],
        "Egfr ligands": ["Egf", "Tgfa", "Hbegf", "Areg", "Btc", "Epgn"],
        "OK-LAM": ["Gpnmb", "Mmp12", "Cd36"],
        "Bad-LAM": ["Trem2", "Ctsd", "Ctsl", "Cd9", "Fabp5"],
        "Tissue-resident VAM": [
            "Lyve1",
            "Cd163",
            "Ednrb",
            "F13a1",
            "Apoe",
            "Cd36",
            "Ccl24",
            "Folr2",
        ],
        "Dividing GH2025": ["Mki67", "Top2a", "Ube2c"],
        "Neural GH2025": ["Siglec1", "Maoa"],
        "Kohda2025": ["Clec4e"],
        "Aging macs GH2025": ["Cxcl13", "C3", "Cd55", "Ly6e", "Ccl8", "Colq", "Mmp9"],
        "Other GH2025": ["Cd209a", "Itgam", "Itgax", "Zbtb46", "Spp1", "Rsad2"],
        "Lipophagy": [
            "Map1lc3b",
            "Sqstm1",
            "Atg16l1",
            "Sra1",
            "Cd36",
            "Sirt6",
            "Pnpla2",
            "Acadm",
            "B4galnt1",
            "Hadh",
            "Inpp5d",
            "Soat1",
            "Ip6k1",
            "Hacl1",
            "Echs1",
            "Dgkz",
            "Hadhb",
        ],
        "Efferocytosis": [
            "Mertk",
            "Lrp1",
            "Havcr1",
            "Timd4",
            "Tyro3",
            "Axl",
            "Cd36",
            "Cx3cr1",
            "Wdfy3",
            "C1qa",
            "C1qb",
            "C1qc",
            "Nr1h3",
            "Nr1h2",
            "Pparg",
        ],
    }  # Ereg not present in filtered sample
    unique_markers = pd.Series(
        list(set([i for j in mac_paper_markers.values() for i in j]))
    )
    assert np.all(unique_markers.isin(adata_macro.var_names))

    mac_markers = {
        "Ccr genes": [
            "Ccr1",
            "Ccr2",
            "Ccr3",
            "Ccr4",
            "Ccr5",
            "Ccr6",
            "Ccr7",
            "Ccr9",
            "Ccr10",
        ],
        "Common IFN response": [
            "Irf9",
            "Irf7",
            "Ifi35",
            "Ifnar2",
            "Isg20",
            "Isg15",
            "Ifit1",
            "Ifit2",
            "Ifih1",
            "Ifnar1",
            "Ifngr2",
            "Cxcl9",
            "Oas1a",
            "Oas1b",
            "Mx1",
        ],
        "Type I IFN": ["Cxcl10", "Isg15", "Mx1", "Irf3", "Irf7", "Il6", "Setdb2"],
        "Type II IFN": ["Irf1", "Stat1", "Mx1", "Bst2", "Mafb"],
        "TGFb response": [
            "Pfkl",
            "Arg1",
            "Cebpa",
            "Id3",
            "Retnla",
            "F13a1",
            "Tgfbr1",
            "Pdcd1",
            "Smad3",
            "Smad7",
            "Cx3cr1",
            "Gcnt2",
            "Pmepa1",
            "Runx3",
            "Axl",
            "F11r",
            "Fcrls",
            "Apoe",
            "Mmp14",
        ],
        "TNF response": [
            "Nfkb1",
            "Nfkbia",
            "Jun",
            "Fos",
            "Mapk1",
            "Mapk3",
            "Mapk8",
            "Mapk8",
            "Mapk10",
            "Srebf2",
            "Il1a",
            "Il1b",
            "Il18",
            "Ccl5",
        ],
        "Type 1 activation": [
            "Cd86",
            "Cd40",
            "Cxcl16",
            "Cxcl9",
            "Il15ra",
            "Il17ra",
            "Icam1",
            "Vcam1",
            "Nfkb1",
            "Rela",
            "Rps6ka2",
            "Tank",
            "Ripk2",
            "Il15",
            "Il23a",
            "Irf1",
            "Slc7a2",
            "Slc12a4",
            "Slc1a4",
            "Slc39a14",
            "Slc3a2",
            "Slc4a7",
            "Csf2ra",
            "Csf2rb",
            "Csf2rb2",
        ],
        "Type 2 activation": [
            "Folr2",
            "Mertk",
            "Arg1",
            "Chil3",
            "Retnla",
            "Aldh1a2",
            "Lpxn",
            "Dhrs3",
            "Mical1",
            "Dnmt3a",
            "Jun",
            "Gab1",
            "P2ry1",
            "Gab2",
            "Glul",
            "Cpt1a",
            "Acadm",
            "B4galnt1",
            "Hadh",
            "Inpp5d",
            "Soat1",
            "Ip6k1",
            "Hacl1",
            "Echs1",
            "Dgkz",
            "Hadhb",
            "Parp1",
            "Ptpn22",
            "Cd300a",
            "Cd84",
            "Csf1",
            "Csf1r",
        ],
        "Non-specific activation": [
            "Ccl7",
            "Ccl17",
            "Ccl22",
            "Ccl24",
            "Cd38",
            "Cd44",
            "Cd83",
            "Csf3r",
        ],
        "Selected markers": ["Gas6", "Il1rn", "S100a4", "Maf1", "Col1a2", "Col5a2"],
    }  # Ucp1, 'Csf3' not present in filtered sample
    unique_markers = pd.Series(list(set([i for j in mac_markers.values() for i in j])))
    assert np.all(unique_markers.isin(adata_macro.var_names))

    fb_markers = {
        "Universal FB markers": [
            "Pdgfra",
            "Dcn",
            "Cxcl12",
            "C3",
            "Col1a1",
            "Fbln5",
            "Lpar1",
            "Ly6a",
        ],
        "Progenitors": ["Cd55", "Pi16", "Dpp4", "Tcf4", "Dlk1", "Islr", "Fn1", "Itga5"],
        "Myofibroblasts": ["Acta2", "Cthrc1"],
        "Adipose-committed": ["Cd24a", "Itgb1", "Icam1", "Apoe", "Fabp5", "Fabp4"],
        "Chondrocyte-committed": ["Sox9", "Col2a1"],
        "Reticular": [
            "Thy1",
            "Wt1",
            "Fap",
            "Acta2",
            "Cnn1",
            "Tgm2",
            "Mgp",
            "Bst1",
            "Mcam",
        ],
        "Angiogenic": ["Cxcl12", "Vegfa", "Vegfd", "Timp1"],
        "Egfr ligands": ["Egf", "Tgfa", "Hbegf", "Areg", "Btc", "Ereg", "Epgn"],
        "Dividing genes": ["Mki67", "Top2a", "Ube2c"],
        "Complement genes": ["C2", "C7", "C4a", "C4b", "Atg7"],
        "Ungrouped": ["Saa3", "Lcn2", "Npnt", "Ces1d", "Cxcl13", "Cxcl12", "Ccl2"],
        "TGFb-resp low": ["Smad2", "Smad3", "Smad4", "Serpine1", "Id1", "Id2", "Id3"],
        "TGFb-resp high": [
            "Fn1",
            "Postn",
            "Cthrc1",
            "Vim",
            "Col1a1",
            "Col1a2",
            "Serpinh1",
        ],
        "IFN genes": [
            "Fasn",
            "Irs1",
            "Insr",
            "Irf3",
            "Irf4",
            "Ifi47",
            "Irf1",
            "Ifitm1",
            "Ifitm2",
            "Ifitm3",
            "Isg15",
            "Cxcl9",
            "Cxcl10",
        ],
        "General IFN/TGFb": [
            "Ifngr1",
            "Ifnar1",
            "Tgfbr2",
            "Il33",
            "Il31ra",
            "Il17ra",
            "Il17rb",
        ],
    }

    unique_markers = pd.Series(list(set([i for j in fb_markers.values() for i in j])))
    assert np.all(unique_markers.isin(adata_fibro.var_names))

    fb_paper_markers = {
        "Emont genes": ["Pde11a", "Aldh1a3", "Mgp", "Tenm3", "Epha3", "Frem1"],
        "Merrick genes": ["Dpp4", "Icam1", "F3", "Wnt6", "Thbs4", "Egfl6"],
        "Burl genes": ["Icam1", "Dpp4", "Lipe", "Top2a"],
        "Hepler genes": ["Fabp4", "Cd36", "Ly6c1"],
        "Sarvari genes": ["Foxp2", "Cd36", "Ebf2", "Klf4"],
        "Wang genes": ["Apoe", "Lifr", "Cd55", "Mfap4"],
        "Kohda genes": ["Osmr"],
    }
    unique_markers = pd.Series(
        list(set([i for j in fb_paper_markers.values() for i in j]))
    )
    assert np.all(unique_markers.isin(adata_fibro.var_names))

    ################################
    ################################
    ################################

    start = stopwatch(setting=2)

    CLUSTER_KEY = "leiden_macro"
    DE_KEY = "macro_DEGs"
    SELECT_RES = [0.6]

    for res in SELECT_RES:
        task_start = stopwatch("Plotting", start)
        res_key = f"{CLUSTER_KEY}_{res:.1f}"
        plot_cluster_violinplot(
            adata_macro,
            "Diet",
            res_key,
            {
                "Arg1 FALC Macs": mac_paper_markers | mac_markers,
            },
        )
        plt.savefig(PLOTDIR / "cluster_violins" / f"{res_key}-Diet.jpg")
        stopwatch(start=task_start, setting=1)

    CLUSTER_KEY = "leiden_fibro"
    DE_KEY = "fibro_DEGs"
    SELECT_RES = [0.5, 0.7]

    for res in SELECT_RES:
        res_key = f"{CLUSTER_KEY}_{res:.1f}"
        plot_cluster_violinplot(
            adata_fibro,
            "Diet",
            res_key,
            fb_markers | fb_paper_markers,
        )
        plt.savefig(PLOTDIR / "cluster_violins" / f"{res_key}-Diet.jpg")

    os.system("zip -r ./plots/cluster_violins.zip ./plots/cluster_violins")
