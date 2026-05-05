# BiocManager::install("POWSC") 
# suppressMessages(library(POWSC))
# library(here)
# obj <- readRDS(here("data", "processed", "single_cell", "Emont2022", "cleaned.rds"))
# 
# ct_counts <- obj$cell_type__custom
# ct_counts <- (ct_counts / sum(ct_counts)) %>% sort() %>% round(3)
# 
# 
# sim_size = 500 # use a large sample is preferrable
# cell_per = c(0.35, 0.25, 0.2)#, 0.08, 0.08, 0.04)
# celltypes = c("macrophage", "ASPC", "adipocyte")#, "mesothelial", "epithelial", "endothelial/mural", "other immune")
# col = obj@meta.data[c("Identifier", "cell_type__custom")]
# sce <- as.SingleCellExperiment(obj)
# sce <- sce[VariableFeatures(obj)]
# exprs = assays(sce)$counts
# 
# estParas_set = NULL
# for (cp in celltypes){
#   print(cp)
#   ix = grep(cp, col$cell_type__custom) %>% sample(2000)
#   tmp_paras = Est2Phase(exprs[, ix] %>% as.matrix())
#   estParas_set[[cp]] = tmp_paras
# }
# 
# 
# data("GSE67835_AB_S7_sce")
# sim_size = 200 # use a large sample is preferrable
# cell_per = c(0.2, 0.3, 0.5)
# col = colData(sce)
# exprs = assays(sce)$counts
# (tb = table(colData(sce)$Patients, colData(sce)$cellTypes))
# # use AB_S7 patient as example and take three cell types: astrocytes hybrid and neurons
# estParas_set = NULL
# celltypes = c("oligodendrocytes", "hybrid", "neurons")
# for (cp in celltypes){
#   print(cp)
#   ix = intersect(grep(cp, col$cellTypes), grep("AB_S7", col$Patients))
#   tmp_mat = exprs[, ix]
#   tmp_paras = Est2Phase(tmp_mat)
#   estParas_set[[cp]] = tmp_paras
# }
# 
# ######### 
# #########  Simulation part
# ######### 
# sim = SimulateMultiSCEs(n = sim_size, estParas_set = estParas_set, multiProb = cell_per)
# 
# ######### 
# #########  DE analysis part
# ######### 
# DE_rslt = NULL
# for (comp in names(sim)){
#   tmp = runDE(sim[[comp]]$sce, DE_Method = "MAST")
#   DE_rslt[[comp]] = tmp
# }
# 
# ######### 
# ######### Summarize the power result
# ######### 
# pow_rslt = pow1 = pow2 = pow1_marg = pow2_marg = NULL
# TD = CD = NULL
# for (comp in names(sim)){
#   tmp1 = Power_Disc(DE_rslt[[comp]], sim[[comp]])
#   tmp2 = Power_Cont(DE_rslt[[comp]], sim[[comp]])
#   TD = c(TD, tmp2$TD); CD = c(CD, tmp2$CD)
#   pow1_marg = c(pow1_marg, tmp1$power.marginal)
#   pow2_marg = c(pow2_marg, tmp2$power.marginal)
#   pow_rslt[[comp]] = list(pow1 = tmp1, pow2 = tmp2)
#   pow1 = rbind(pow1, tmp1$power)
#   pow2 = rbind(pow2, tmp2$power)
# }
# 
# ######### 
# ######### Demonstrate the result by heatmap
# ######### 
# library(RColorBrewer); library(pheatmap)
# breaksList = seq(0, 1, by = 0.01)
# colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList))
# dimnames(pow1) = list(names(sim), names(tmp1$CD))
# dimnames(pow2) = list(names(sim), names(tmp2$CD))
# pheatmap(pow2, display_numbers = TRUE, color=colors, show_rownames = TRUE,
#          cellwidth = 30, cellheight = 40, legend = TRUE,
#          border_color = "grey96", na_col = "grey",
#          cluster_row = FALSE, cluster_cols = FALSE,
#          breaks = seq(0, 1, 0.01),
#          main = "")



devtools::install_github("immunogenomics/scpost")
library(scpost)
library(Seurat)
library(here)

# using example Emont2022 prior, we estimate 20k cell/sample from 60k cell/lane (3 lanes/sample from 15samples over 5 lanes)

obj <- readRDS(here("data", "processed", "single_cell", "Emont2022", "cleaned.rds"))

raFib_freqEstimates <- estimateFreqVar(
  meta = obj@meta.data, clusCol = 'cell_type__custom', sampleCol = 'Identifier', logCov = TRUE)

raFib_pcEstimates <- estimatePCVar(
  pca = Embeddings(object = obj, reduction = "PCA_hvg"), npcs = 20,
  meta = obj@meta.data, clusCol = 'cell_type__custom', sampleCol = 'Identifier', batchCol = 'donor_id')

centroids <- raFib_pcEstimates$centroids
pc_cov_list <- raFib_pcEstimates$pc_cov_list
batch_vars <- raFib_pcEstimates$batch_vars
sample_vars <- raFib_pcEstimates$sample_vars

set.seed(123)

# Set the number of samples, number of cells per sample, and create batch structure
ncases <- 5
nctrls <- 5
nbatches <- 2
batchStructure <- distribSamples(ncases = ncases, nctrls = nctrls, nbatches = nbatches)
ncells <- rep(5000, times = ncases + nctrls)
names(ncells) <- batchStructure$sample_names

params <- createParamTable(
  nreps = 10,
  clus = "clusASPC",
  fc = 5,
  ncases = ncases,
  nctrls = nctrls,
  nbatches = nbatches,
  b_scale = 1,
  s_scale = 1,
  cf_scale = 1,
  res_use = 0.6,
  cond_induce = "cases",
  save_path = file.path(getwd(), "power_analysis")
)

suppressWarnings({
  lapply(seq(nrow(params)), function(x){
    simDataset.withMASC(
      save_path = params[x, 'save_path'],
      rep = params[x, 'rep'],
      seed = params[x, 'seed'],
      ncases = params[x, 'ncases'],
      nctrls = params[x, 'nctrls'],
      nbatches = params[x, 'nbatches'],
      batchStructure = batchStructure,
      ncells = ncells,
      centroids = centroids,
      pc_cov_list = pc_cov_list,
      batch_vars = batch_vars,
      b_scale = params[x, 'b_scale'],
      sample_vars = sample_vars,
      s_scale = params[x, 's_scale'],
      cfcov = raFib_freqEstimates$cfcov,
      cf_scale = params[x, 'cf_scale'],
      meanFreqs = raFib_freqEstimates$meanFreq,
      clus = params[x, 'clus'],
      fc = params[x, 'fc'],
      cond_induce = params[x, 'cond_induce'],
      res_use = params[x, 'res_use'], 
      mc.cores = 1,
      clusterData = TRUE,
      returnPCs = FALSE
    )
  })
})

dir <- file.path(getwd(), "power_analysis")
filenames <- list.files(path = dir,
                        full.names = T,
                        pattern = '*res') %>% basename
resTables <- lapply(filenames, function(x){
  readRDS(file.path(dir, x))[["res"]]
})

getPowerFromRes(
  resFiles = filenames,
  resTables = resTables,
  threshold = 0.05,
  z = 1.96,
  stratByClus = FALSE
)
