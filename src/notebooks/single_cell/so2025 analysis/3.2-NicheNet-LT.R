# general
options(max.print = 300)
options(timeout = 300)
options(stringsAsFactors = FALSE)
library(dplyr)
library(here)
library(Seurat)
library(nichenetr)
source(here("src", "single_cell", "nichenet.R"))

clean_names <- function(cell.name) {
  cell.name <- gsub("\\?\\?\\?", "unknown", cell.name)
  cell.name <- gsub("\\+\\/", "+-", cell.name)
  cell.name <- gsub("\\/", "+", cell.name)
  return(cell.name)
}

# Prepare inputs
print("Read in seurat object...")
filename <- "RDS_annotated_eWAT.rds"
version.name <- "-eWAT"
system.time(seurat_object <- readRDS(here("data","processed", "single cell", "5_analysis", filename)))
celltype_col <- "cell_type_hires"

# Preprocess inputs
print("Start preparing inputs...")
nichenet.inputs <- prepNicheNet(seurat_object, celltype_col)

print("Pulling celltypes of interest...")
# celltypes of interest
celltypes <- c(nichenet.inputs$cell.types[grepl("fibro", nichenet.inputs$cell.types)],
               nichenet.inputs$cell.types[grepl("adipo", nichenet.inputs$cell.types)],
               nichenet.inputs$cell.types[grepl("macro", nichenet.inputs$cell.types)])
celltypes <- celltypes[!grepl("SAT", celltypes)]
# nichenet_results <- list(celltypes=celltypes)

# run NicheNet & plots
for (celltype in celltypes) {
  print("#############")
  print(paste("Beginning cell type:",celltype))

  nichenet_results <- list()

  print("Running NicheNet...")
  nichenet_results$nichenet.outs <- runNicheNet(celltype, nichenet.inputs)
  
  print("Running Agnostic Plots...")
  system.time(nichenet_results$agnostic.plots <- plotNicheNet_Agnostic(
    nichenet.inputs,
    nichenet_results$nichenet.outs,
    n_ligands = 40))

  print("Running Sender-Specific Plots...")
  system.time(nichenet_results$sender.plots <- plotNicheNet_SenderSpecific(
    nichenet.inputs,
    nichenet_results$nichenet.outs,
    n_ligands = 20))
  
  print("Beginning to save plots!")
  savedir <- here("data", "processed", "single cell", "5_analysis", paste0("nichenet",version.name), clean_names(celltype))
  dir.create(savedir, showWarnings = F, recursive=T)
  
  # agnostic plots
  ggsave(
    paste0(savedir, "/agnostic-","ligand_qc.jpg"),
    plot=nichenet_results$agnostic.plots[[1]],
    width=5, height=3)
  ggsave(
    paste0(savedir, "/agnostic-","results.jpg"),
    plot=nichenet_results$agnostic.plots[[2]],
    width=50, height=10, limitsize = FALSE)
  
  # sender-specific plots
  ggsave(
    paste0(savedir, "/sender-","ligand_qc.jpg"),
    plot=nichenet_results$sender.plots[[1]],
    width=5, height=3)
  ggsave(
    paste0(savedir, "/sender-","results.jpg"),
    plot=nichenet_results$sender.plots[[2]],
    width=50, height=10, limitsize = FALSE)
  # ggsave(
  #   paste0(savedir, "/comparison.jpg"),
  #   plot=nichenet_results$sender.plots[[3]],
  #   width=10, height=10)
  
  print("#############")
}

# saveRDS(nichenet_results, here("data", "processed", "single cell", "5_analysis", "NicheNet_hires_results.rds"))