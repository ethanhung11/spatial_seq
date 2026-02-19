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
# TODO change name & filename
print("Read in seurat object...")
filename <- "RDS_object.RDS"
celltype_col <- "celltype"
descriptor <- "eWAT"
input_file <- here("data", "processed", "single_cell", "combined", filename)
save_dir <- here("data", "processed", "single_cell", "combined", "tools", "nichenet", descriptor)
RUN <- TRUE

# Preprocess inputs
print("Start preparing inputs...")
system.time(seurat_object <- readRDS(input_file))
if (RUN == TRUE) {
  nichenet.inputs <- prepNicheNet(seurat_object, celltype_col)
  saveRDS(nichenet.inputs, here(save_dir, "nichenet_inputs.RDS"))
} else {
  nichenet.inputs <- readRDS(nichenet.results, here(save_dir, "nichenet_inputs.RDS"))
}


print("Pulling celltypes of interest...")
# celltypes of interest
celltypes <- c(nichenet.inputs$cell.types[grepl("fibro", nichenet.inputs$cell.types)],
               nichenet.inputs$cell.types[grepl("adipo", nichenet.inputs$cell.types)],
               nichenet.inputs$cell.types[grepl("macro", nichenet.inputs$cell.types)])

# run NicheNet & plots
for (celltype in celltypes) {
  print("#############")
  print(paste("Beginning cell type:",celltype))
  celltype_save_dir <- here(save_dir,celltype)
  dir.create(celltype_save_dir, showWarnings = F, recursive=T)

  nichenet.results <- list()
  if (RUN == TRUE) {
    print("Running NicheNet...")
    nichenet.results$nichenet.outs <- runNicheNet(celltype, nichenet.inputs)
    
    print("Running Agnostic Plots...")
    system.time(nichenet.results$agnostic.plots <- plotNicheNet_Agnostic(
      nichenet.inputs,
      nichenet.results$nichenet.outs,
      n_ligands = 40))

    print("Running Sender-Specific Plots...")
    system.time(nichenet.results$sender.plots <- plotNicheNet_SenderSpecific(
      nichenet.inputs,
      nichenet.results$nichenet.outs,
      n_ligands = 20))

    saveRDS(nichenet.results, here(celltype_save_dir, paste0("nichenet_",celltype,".RDS")))
  } else {
    nichenet.results <- readRDS(nichenet.results, paste0("nichenet_",celltype,".RDS"))
  }


  print("Beginning to save plots!")
  # agnostic plots
  ggsave(
    paste0(celltype_save_dir, "/agnostic-","ligand_qc.jpg"),
    plot=nichenet.results$agnostic.plots[[1]],
    width=5, height=3)
  ggsave(
    paste0(celltype_save_dir, "/agnostic-","results.jpg"),
    plot=nichenet.results$agnostic.plots[[2]],
    width=50, height=10, limitsize = FALSE)
  
  # sender-specific plots
  ggsave(
    paste0(celltype_save_dir, "/sender-","ligand_qc.jpg"),
    plot=nichenet.results$sender.plots[[1]],
    width=5, height=3)
  ggsave(
    paste0(celltype_save_dir, "/sender-","results.jpg"),
    plot=nichenet.results$sender.plots[[2]],
    width=50, height=10, limitsize = FALSE)
  ggsave(
    paste0(celltype_save_dir, "/comparison.jpg"),
    plot=nichenet.results$sender.plots[[3]],
    width=10, height=10)
  
  print("#############")
}
