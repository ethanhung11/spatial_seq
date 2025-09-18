# general
options(max.print = 300)
options(stringsAsFactors = FALSE)
library(tidyverse)
library(Seurat)
library(nichenetr)

prepNicheNet <- function(seurat_object, celltype_col, condition_oi="HFD", condition_reference="LFD", organism="mouse") {
  # Get DBs
  print("Getting databases...")
  if(organism == "human"){
    lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
    ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
    weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))
  } else if(organism == "mouse"){
    lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
    ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
    weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
  } else {
    stop(paste0("organism '",organism,"' is not allowed. Must use 'human' or 'mouse'"))
  }
  lr_network <- lr_network %>% distinct(from, to)
  
  # # View DBs
  # head(lr_network)
  # ligand_target_matrix[1:5,1:5]
  # head(weighted_networks$lr_sig)
  # head(weighted_networks$gr)
  
  # get cell types & assign them
  print("Getting available celltypes and processing Seurat obj...")
  labels <- seurat_object@meta.data[celltype_col][,] %>% as.factor()
  Idents(seurat_object) <- labels
  names(labels) <- Cells(seurat_object)
  cell.types <- labels %>% levels()
  
  return(list(
    seurat_object=seurat_object,
    celltype_col=celltype_col,
    cell.types=cell.types,
    labels=labels,
    condition_oi=condition_oi,
    condition_reference=condition_reference,
    lr_network=lr_network,
    ligand_target_matrix=ligand_target_matrix,
    weighted_networks=weighted_networks))
}


runNicheNet <- function(receiver, nichenet.ins) {
  # ACCESS NICHENET INS
  seurat_object <- nichenet.ins$seurat_object
  celltype_col <- nichenet.ins$celltype_col
  cell.types <- nichenet.ins$cell.types
  labels <- nichenet.ins$labels
  condition_oi <- nichenet.ins$condition_oi
  condition_reference <- nichenet.ins$condition_reference
  lr_network <- nichenet.ins$lr_network
  ligand_target_matrix <- nichenet.ins$ligand_target_matrix
  weighted_networks <- nichenet.ins$weighted_networks
  
  # Get the expressed genes of the receiver cell type separately here
  expressed_genes_receiver <- get_expressed_genes(receiver, GetAssayData(seurat_object), labels, 0.05)
  all_receptors <- unique(lr_network$to)  
  expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)
  
  # Get the expressed genes of every sender cell type separately here
  senders <- cell.types[cell.types!=receiver]
  potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()
  list_expressed_genes_sender <- senders %>% lapply(get_expressed_genes, GetAssayData(seurat_object), labels, 0.05)
  expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()
  sender_potential_ligands <- intersect(potential_ligands, expressed_genes_sender) 
  
  # Find genes of interest. Can also take from AnnData object.
  # calculate DEGs
  receiver_idx <- names(labels)[labels == celltype]
  seurat_obj_receiver <- subset(seurat_object, cells = receiver_idx)
  DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,
                                    ident.1 = condition_oi, ident.2 = condition_reference,
                                    group.by = "Condition",
                                    min.pct = 0.05) %>% rownames_to_column("gene")
  
  # Define geneset using DEGs
  geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
  geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
  # Define background
  background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
  
  # Run NicheNet
  ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                                 background_expressed_genes = background_expressed_genes,
                                                 ligand_target_matrix = ligand_target_matrix,
                                                 potential_ligands = potential_ligands)
  ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = row_number())
  
  return(list(
    receiver=receiver,
    senders=senders,
    geneset_oi=geneset_oi,
    ligand_activities=ligand_activities,
    expressed_receptors=expressed_receptors,
    sender_potential_ligands=sender_potential_ligands))
}


plotNicheNet_Agnostic <- function(nichenet.ins, nichenet.outs, n_ligands) {
  # ACCESS NICHENET INS
  lr_network <- nichenet.ins$lr_network
  ligand_target_matrix <- nichenet.ins$ligand_target_matrix
  weighted_networks <- nichenet.ins$weighted_networks
  
  # ACCESS NICHENET OUTS
  receiver <- nichenet.outs$receiver
  ligand_activities <- nichenet.outs$ligand_activities
  geneset_oi <- nichenet.outs$geneset_oi
  expressed_receptors <- nichenet.outs$expressed_receptors
  
  
  # LIGAND CUTOFFS
  p_hist_lig_activity <- ggplot(ligand_activities, aes(x=aupr_corrected)) + 
    geom_histogram(color="black", fill="darkorange", bins = 100)  + 
    geom_vline(aes(xintercept=min(ligand_activities %>% top_n(n_ligands, aupr_corrected) %>% pull(aupr_corrected))),
               color="red", linetype="dashed", size=1) + 
    labs(x="ligand activity (PCC)", y = "# ligands") +
    theme_classic()
  
  
  # TOP LIGANDS
  best_upstream_ligands <- ligand_activities %>%
    top_n(n_ligands, aupr_corrected) %>%
    arrange(-aupr_corrected) %>% pull(test_ligand)
  vis_ligand_aupr <- ligand_activities %>%
    filter(test_ligand %in% best_upstream_ligands) %>%
    column_to_rownames("test_ligand") %>%
    select(aupr_corrected) %>% arrange(aupr_corrected) %>%
    as.matrix(ncol = 1)
  
  p_ligand_aupr <- make_heatmap_ggplot(vis_ligand_aupr,
                                       "Prioritized ligands", "Ligand activity", 
                                       legend_title = "AUPR", color = "darkorange") + 
    theme(axis.text.x.top = element_blank())
  
  
  # TOP TARGETS
  active_ligand_target_links_df <- best_upstream_ligands %>%
    lapply(get_weighted_ligand_target_links,
           geneset = geneset_oi,
           ligand_target_matrix = ligand_target_matrix,
           n = 50) %>%
    bind_rows() %>% drop_na()
  
  print(paste("ligand-target links:",nrow(active_ligand_target_links_df)))
  active_ligand_target_links <- prepare_ligand_target_visualization(
    ligand_target_df = active_ligand_target_links_df,
    ligand_target_matrix = ligand_target_matrix,
    cutoff = 0.1) 
  
  print(paste("Resulting links after cutoff:",nrow(active_ligand_target_links)))
  order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
  order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))
  
  vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])
  
  p_ligand_target <- make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                                         color = "purple", legend_title = "Regulatory potential") +
    scale_fill_gradient2(low = "whitesmoke",  high = "purple")
  
  
  # TOP RECEPTORS
  ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
    best_upstream_ligands, expressed_receptors,
    lr_network, weighted_networks$lr_sig) 
  
  vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
    ligand_receptor_links_df,
    best_upstream_ligands,
    order_hclust = "receptors") 
  
  p_ligand_receptor <- make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                                           y_name = "Ligands", x_name = "Receptors",  
                                           color = "mediumvioletred", legend_title = "Prior interaction potential")
  
  # PLOT TOGETHER
  par(mar=c(0.1,0.1,0.1,0.1))
  ratio <- c(ncol(vis_ligand_aupr)+6,
             nrow(vis_ligand_receptor_network),
             ncol(vis_ligand_target))
  
  figures_without_legend <- cowplot::plot_grid(
    p_ligand_aupr + theme(legend.position = "none"),
    p_ligand_receptor + theme(legend.position = "none",
                              axis.title.y = element_blank()),
    p_ligand_target + theme(legend.position = "none",
                            axis.title.y = element_blank()),
    align = "hv",
    nrow = 1,
    rel_widths = ratio)
  
  legends <- cowplot::plot_grid(
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_receptor)),
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target)),
    nrow = 1,
    align = "h", scale=0.5, rel_widths = ratio)
  
  combined_plot <-  cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
  
  return(list(ligand_qc=plot(p_hist_lig_activity),result=plot(combined_plot)))
}


plotNicheNet_SenderSpecific <- function(nichenet.ins, nichenet.outs, n_ligands) {
  # ACCESS NICHENET INS
  seurat_object <- nichenet.ins$seurat_object
  celltype_col <- nichenet.ins$celltype_col
  cell.types <- nichenet.ins$cell.types
  labels <- nichenet.ins$labels
  condition_oi <- nichenet.ins$condition_oi
  condition_reference <- nichenet.ins$condition_reference
  lr_network <- nichenet.ins$lr_network
  ligand_target_matrix <- nichenet.ins$ligand_target_matrix
  weighted_networks <- nichenet.ins$weighted_networks
  
  # ACCESS NICHENET OUTS
  receiver <- nichenet.outs$receiver
  senders <- nichenet.outs$senders
  geneset_oi <- nichenet.outs$geneset_oi
  expressed_receptors <- nichenet.outs$expressed_receptors
  ligand_activities <- nichenet.outs$ligand_activities
  sender_potential_ligands <- nichenet.outs$sender_potential_ligands
  
  # LIGAND CUTOFFS
  ligand_activities_all <- ligand_activities
  ligand_activities <- ligand_activities %>% filter(test_ligand %in% sender_potential_ligands)
  best_upstream_ligands <- ligand_activities %>% top_n(n_ligands, aupr_corrected) %>% arrange(-aupr_corrected) %>%
    pull(test_ligand) %>% unique()
  
  p_hist_lig_activity <- ggplot(ligand_activities, aes(x=aupr_corrected)) + 
    geom_histogram(color="black", fill="darkorange", bins = 100)  + 
    geom_vline(aes(xintercept=min(ligand_activities %>% top_n(n_ligands, aupr_corrected) %>% pull(aupr_corrected))),
               color="red", linetype="dashed", size=1) + 
    labs(x="ligand activity (PCC)", y = "# ligands") +
    theme_classic()
  
  # TOP LIGANDS
  ligand_aupr_matrix <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
    column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected)
  vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1) 
  
  p_ligand_aupr <- make_heatmap_ggplot(
    vis_ligand_aupr,
    "Prioritized ligands", "Ligand activity", 
    legend_title = "AUPR", color = "darkorange") + 
    theme(axis.text.x.top = element_blank())
  
  # TOP TARGETS
  active_ligand_target_links_df <- best_upstream_ligands %>%
    lapply(get_weighted_ligand_target_links,
           geneset = geneset_oi,
           ligand_target_matrix = ligand_target_matrix,
           n = 50) %>%
    bind_rows() %>% drop_na()
  
  active_ligand_target_links <- prepare_ligand_target_visualization(
    ligand_target_df = active_ligand_target_links_df,
    ligand_target_matrix = ligand_target_matrix,
    cutoff = 0.33) 
  
  order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
  order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))
  
  vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])
  
  p_ligand_target <- make_heatmap_ggplot(
    vis_ligand_target,
    "Prioritized ligands",
    "Predicted target genes",
    color = "purple",
    legend_title = "Regulatory potential") +
    scale_fill_gradient2(low = "whitesmoke",  high = "purple")
  
  # TOP RECEPTORS
  ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
    best_upstream_ligands, expressed_receptors,
    lr_network, weighted_networks$lr_sig) 
  
  vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
    ligand_receptor_links_df,
    best_upstream_ligands,
    order_hclust = "receptors") 
  
  p_ligand_receptor <- make_heatmap_ggplot(
    t(vis_ligand_receptor_network), 
    y_name = "Ligands", x_name = "Receptors",  
    color = "mediumvioletred", legend_title = "Prior interaction potential")
  
  # SENDER-SPECIFIC INFO
  p_dotplot <- DotPlot(
    subset(seurat_object, cells = names(labels)[labels %in% senders]),
    features = rev(best_upstream_ligands), cols = "RdYlBu") +
    coord_flip() +
    scale_y_discrete(position = "right")
  
  celltype_order <- levels(Idents(seurat_object)) 
  DE_table_top_ligands <- lapply(
    celltype_order[celltype_order %in% senders],
    get_lfc_celltype, 
    seurat_obj = seurat_object,
    condition_colname = "Condition",
    condition_oi = condition_oi,
    condition_reference = condition_reference,
    celltype_col = celltype_col,
    min.pct = 0, logfc.threshold = 0,
    features = best_upstream_ligands 
  ) 
  
  DE_table_top_ligands <- DE_table_top_ligands %>%  reduce(., full_join) %>% 
    column_to_rownames("gene")
  vis_ligand_lfc <- as.matrix(DE_table_top_ligands[rev(best_upstream_ligands), , drop = FALSE])
  p_lfc <- make_threecolor_heatmap_ggplot(
    vis_ligand_lfc,
    "Prioritized ligands", "LFC in Sender",
    low_color = "midnightblue", mid_color = "white",
    mid = median(vis_ligand_lfc), high_color = "red",
    legend_title = "LFC")
  
  # COMPARISON
  comparison_plot <- (make_line_plot(ligand_activities = ligand_activities_all,
                                     potential_ligands = ligand_activities$test_ligand,
                                     ranking_range=c(1,40)) +
                        theme(plot.title = element_text(size=11, hjust=0.1, margin=margin(0, 0, -5, 0))))
  
  # COMBINED PLOT
  par(mar=c(0.1,0.1,0.1,0.1))
  ratio <- c(ncol(vis_ligand_aupr)+6,
             ncol(vis_ligand_lfc)+7,
             ncol(vis_ligand_lfc)+8,
             nrow(vis_ligand_receptor_network),
             ncol(vis_ligand_target))
  
  figures_without_legend <- cowplot::plot_grid(
    p_ligand_aupr + theme(legend.position = "none"),
    p_dotplot + theme(legend.position = "none",
                      axis.ticks = element_blank(),
                      axis.title.y = element_blank(),
                      axis.title.x = element_text(size = 12),
                      axis.text.y = element_text(size = 9),
                      axis.text.x = element_text(size = 9,  angle = 90, hjust = 0)) +
      ylab("Expression in Sender"),
    p_lfc + theme(legend.position = "none",
                  axis.title.y = element_blank()),
    p_ligand_receptor + theme(legend.position = "none",
                              axis.title.y = element_blank()),
    p_ligand_target + theme(legend.position = "none",
                            axis.title.y = element_blank()),
    align = "hv", nrow = 1, rel_widths = ratio)
  
  
  legends <- cowplot::plot_grid(
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
    ggpubr::as_ggplot(ggpubr::get_legend(p_dotplot)),
    ggpubr::as_ggplot(ggpubr::get_legend(p_lfc)),
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_receptor)),
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target)),
    nrow = 1, align = "h", scale=0.5, rel_widths = ratio)
  
  combined_plot <-  cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
  
  return(list(ligand_qc=plot(p_hist_lig_activity), result=plot(combined_plot), comparison=plot(comparison_plot)))
}