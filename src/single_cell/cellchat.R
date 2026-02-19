library(CellChat)
library(NeuronChat)
library(Seurat)

prepCellChat <- function(
    seurat.obj,
    celltype_col="celltype",
    sample_col=NULL,
    db=NULL,
    spatial=FALSE,
    spatial_key="spatial",
    spatial_factors=NULL) {

  Idents(seurat.obj) <- unlist(seurat.obj[[celltype_col]])
  
  if (!is.null(sample_col)) {
    seurat.obj[["samples"]] <- seurat.obj[[sample_col]]
  }
  
  if (spatial == TRUE) {
    cellchat <- createCellChat(
      object = seurat.obj[["RNA"]]@data,
      meta = seurat.obj@meta.data,
      group.by = celltype_col,
      datatype = "spatial",
      coordinates = Embeddings(seurat.obj,spatial_key),
      spatial.factors = spatial_factors
    )
  } else {
    cellchat <- createCellChat(
      object = seurat.obj[["RNA"]]@data,
      meta = seurat.obj@meta.data,
      group.by = celltype_cols
    )
  }
  
  if (is.null(db)) {
    CellChatDB <- CellChatDB.mouse
    CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")
    cellchat@DB <- CellChatDB.use
  } else {
    cellchat@DB <- db
  }
  
  return(cellchat)
}

runCellChat <- function(cellchat, spatial=FALSE, pathway=FALSE) {
  options(future.globals.maxSize = 1000e6)
  future::plan("multisession", workers = 1)
  print("filtering overexpressed genes -> LR pairs")
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # calculate cell communication
  print("running LR modeling")
  if (spatial == TRUE) {
    cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, 
                                  distance.use = FALSE, interaction.range = 250, scale.distance = NULL,
                                  contact.dependent = TRUE, contact.range = 100)
  } else {
    cellchat <- computeCommunProb(cellchat, type = "triMean")
  }
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  print("computing network centrality measures")
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "net")
  
  # calculate pathway
  if (pathway == TRUE) {
    print("running pathway modeling as well!")
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  }
  
  return(cellchat)
}

mergeCellChat.custom <- function(cellchat.list, pathway=FALSE) {
  cellchat.MERGED <- mergeCellChat(cellchat.list, add.names = names(cellchat.list))
  available_slots <- if (pathway) c("net", "netP") else c("net")
  for (slot in available_slots) {
    cellchat.MERGED <- computeNetSimilarityPairwise(cellchat.MERGED, type = "functional", slot.name = slot)
    cellchat.MERGED <- netEmbedding(cellchat.MERGED, type = "functional", slot.name = slot)
    cellchat.MERGED <- netClustering(cellchat.MERGED, type = "functional", slot.name = slot)
  }
  return(cellchat.MERGED)
}

customCellChatDB <- function(lr_data, pathway_name = "Custom", 
                                      annotation = "Secreted Signaling",
                                      evidence = "Custom Database") {
  # Create interaction table
  interaction_df <- data.frame(
    interaction_name = paste0(lr_data$ligand, "_", lr_data$receptor),
    pathway_name = pathway_name,
    ligand = lr_data$ligand,
    receptor = lr_data$receptor,
    agonist = "",
    antagonist = "",
    co_A_receptor = "",
    co_I_receptor = "",
    evidence = evidence,
    annotation = annotation,
    interaction_name_2 = paste(lr_data$ligand, "-", lr_data$receptor),
    stringsAsFactors = FALSE
  )
  
  # Create gene info table
  all_genes <- unique(c(lr_data$ligand, lr_data$receptor))
  gene_info <- data.frame(Symbol = all_genes, stringsAsFactors = FALSE)
  
  # Create empty complex and cofactor tables
  complex_df <- data.frame(
    complex_name = character(0),
    subunit_1 = character(0),
    subunit_2 = character(0),
    subunit_3 = character(0),
    subunit_4 = character(0),
    stringsAsFactors = FALSE
  )
  
  cofactor_df <- data.frame(
    cofactor_name = character(0),
    cofactor_1 = character(0),
    cofactor_2 = character(0),
    cofactor_3 = character(0),
    stringsAsFactors = FALSE
  )
  
  # Create CellChatDB object
  CellChatDB <- list(
    interaction = interaction_df,
    complex = complex_df,
    cofactor = cofactor_df,
    geneInfo = gene_info
  )
  
  return(CellChatDB)
}


runNeuronChat <- function(
    seurat.obj,
    celltype_col="celltype",
    spatial_key=NULL) {
  
  if (is.null(spatial_key)) {
    meta <- seurat.obj@meta.data
  } else {
    meta <- cbind(seurat.obj@meta.data, Embeddings(seurat.obj,spatial_key))
  }
  
  nc_obj <- NeuronChat::createNeuronChat(seurat_object[["RNA"]]@data,
                                         DB = "mouse",
                                         meta = meta,
                                         group.by = unlist(seurat.obj[[celltype_col]]))
  
  nc_obj <- NeuronChat::run_NeuronChat(nc_obj, M = 100)
  return(nc_obj)
}