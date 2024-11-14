library(Seurat)
library(DoubletFinder)
library(clustree)
library(parallel)

options(mc.cores = 1)



################################################################################
### QC & doublet
setwd("/xdisk/mliang1/qqiu/project/others/test/cluster/")
outfile <- "/xdisk/mliang1/qqiu/project/others/test/cluster/huves1SN.cluster.rds"

seurat_data <- Read10X(data.dir = "/xdisk/mliang1/qqiu/project/others/test/cellranger/huves1SN/outs/filtered_feature_bc_matrix/")
seurat_object <- CreateSeuratObject(counts = seurat_data, project = "huves1SN", min.cells = 3, min.features = 200)

seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

### QC - doubletfinder
seurat_object <- NormalizeData(seurat_object)
seurat_object <- FindVariableFeatures(seurat_object)
seurat_object <- ScaleData(seurat_object)
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
seurat_object <- RunUMAP(seurat_object, dims=1:40)

stdv <- seurat_object[["pca"]]@stdev
sum_stdv <- sum(seurat_object[["pca"]]@stdev)
percent_stdv <- (stdv / sum_stdv) * 100
cumulative <- cumsum(percent_stdv)
co1 <- which(cumulative > 90 & percent_stdv < 5)[1]
co2 <- sort(which((percent_stdv[1:length(percent_stdv) - 1] -
                    percent_stdv[2:length(percent_stdv)]) > 0.1),
           decreasing = T)[1] + 1
pc <- min(co1, co2)
pc_list <- c(pc, 20, 30, 40)

para_list <- data.frame(sample_ID = character(),
                       n_cell = numeric(),
                       n_doublet = numeric(),
                       doublet_rate = numeric(),
                       PC = numeric(),
                       pK = numeric())

## try different pcs
for(pci in 1:length(pc_list)){
  
  pc <- pc_list[pci]
  ## pK identification
  sweep.list <- paramSweep_v3(seurat_object, PCs = 1:pc, num.cores = 1)
  sweep.stats <- summarizeSweep(sweep.list)
  bcmvn <- find.pK(sweep.stats)
  
  pK <- as.numeric(as.character(bcmvn$pK))
  BCmetric <- bcmvn$BCmetric
  pK_choose <- pK[which(BCmetric %in% max(BCmetric))]
  
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  
  doublet_rate <- ncol(seurat_object)*8*1e-6
  if((doublet_rate > 0.3) & (doublet_rate > optimal.pk)){
    doublet_rate=optimal.pk
  }
  
  ## Homotypic doublet proportion estimate
  annotations <- seurat_object@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  nExp.poi <- round(doublet_rate * nrow(seurat_object@meta.data)) ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  
  # run DoubletFinder
  seurat_object <- doubletFinder_v3(seu = seurat_object, 
                                   PCs = 1:pc, 
                                   pK = optimal.pk,
                                   nExp = nExp.poi.adj)
  
  metadata <- seurat_object@meta.data
  if(pci==1){
    colnames(metadata)[ncol(metadata)] <- "doublet_pc.optm"
  }else{
    colnames(metadata)[ncol(metadata)] <- paste0("doublet_pc.", pci)
  }
  
  seurat_object@meta.data <- metadata
  
  seurat_object@meta.data <- seurat_object@meta.data[, -which(grepl("pANN", colnames(seurat_object@meta.data)))]
  
  
  para_list[nrow(para_list)+1, ] <- c("huves1SN", nrow(metadata),
                                     table(metadata[,ncol(metadata)])['Doublet'],
                                     doublet_rate,
                                     pc, optimal.pk)
}





################################################################################
### cluster & annotation
black_list <- c(rownames(seurat_object)[grep("^RP[SL][[:digit:]]", rownames(seurat_object))],
               rownames(seurat_object)[grep("^IG", rownames(seurat_object))])
rm_black_list <- setdiff(rownames(seurat_object), black_list)
seurat_object <- subset(seurat_object, nFeature_RNA>200 & 
                         percent.mt<20 & 
                         doublet_pc.2=="Singlet",
                         features=rm_black_list)

### cell cycle score & to keep the signals separating non-cycling and cycling cells; reference: https://satijalab.org/seurat/articles/cell_cycle_vignette.html
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
seurat_object <- CellCycleScoring(seurat_object, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
seurat_object$CC.Difference <- seurat_object$S.Score - seurat_object$G2M.Score

### basic normalization
seurat_object <- NormalizeData(seurat_object)
seurat_object <- FindVariableFeatures(seurat_object)
seurat_object <- ScaleData(seurat_object, features = rownames(seurat_object))
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
ElbowPlot(seurat_object)

nPC <- 20; reso <- c(0.5, 1, 1.5, 2)
seurat_object <- RunUMAP(seurat_object, dims=1:nPC)
seurat_object <- FindNeighbors(seurat_object, dims=1:nPC)
seurat_object <- FindClusters(seurat_object, resolution = reso)
seurat_meta <- seurat_object@meta.data
clustree(seurat_meta, prefix = "RNA_snn_res.")

for(i in c(1)){
  cluster = paste0("RNA_snn_res.", i)
  p = DimPlot(seurat_object, label = T, reduction = "umap", group = cluster)
  print(p)
  # outfile = paste0("rat.ss.MSA.reso_", i)
  # DE_analysis(seurat_object, cluster = cluster, outfile = paste0("../DEG/", outfile))
}


FeaturePlot(seurat_object, marker_list)

umap_dotplot(seurat_object, marker_list, cluster="RNA_snn_res.1", reduction = "umap")

markers <- FindAllMarkers(seurat_object, group.by="RNA_snn_res.1", assay = "RNA", only.pos = TRUE, min.pct = 0.25)
markers$pct.diff = markers$pct.1 - markers$pct.2
write.table(markers,file=gsub("cluster.rds", "allmarker.0.25.long.txt", outfile), sep="\t")

marker_tbl <- markers[markers$p_val_adj<0.05 & markers$avg_log2FC>0.1,] %>% group_by(cluster) %>%
  dplyr::mutate(id = row_number()) %>% as.data.frame() %>%
  select(id, cluster, gene) %>% reshape(., idvar = "id", timevar = "cluster",
                                        v.names="gene", sep=" ", direction = "wide")
colnames(marker_tbl) <- sub("gene", "Cluster", colnames(marker_tbl))
write.table(marker_tbl,file=gsub("cluster.rds", "allmarker.0.25.wide.txt", outfile), sep="\t", row.names = F)

cell_type_annotation <- c("0"="Pericyte", "1"="EC", "2"="Pericyte", "3"="EC", "4"="Fibroblast",
                         "5"="VSMC", "6"="LEC", "7"="Immune")
seurat_object$cell_type <- cell_type_annotation[as.character(seurat_object$RNA_snn_res.1)]
seurat_object$cell_type <- factor(seurat_object$cell_type, levels = c("EC", "LEC", "VSMC", "Pericyte", "Fibroblast", "Macrophage"))

saveRDS(seurat_object, "/xdisk/mliang1/qqiu/project/others/test/cluster/huves1SN.cluster.rds")




################################################################################
### NG-2022 thoracic aorta data process
expr_mtx <- ReadMtx(mtx = "/xdisk/mliang1/qqiu/project/others/test/data/NG-2022/ascending_descending_human_aorta_v1.mtx",
                 cells = "/xdisk/mliang1/qqiu/project/others/test/data/NG-2022/ascending_descending_human_aorta_v1_barcodes.tsv",
                 features = "/xdisk/mliang1/qqiu/project/others/test/data/NG-2022/ascending_descending_human_aorta_v1_genes.tsv")

metadata <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/NG-2022/ascending_descending_human_aorta_metadata.txt", 
                       header = TRUE, sep = "\t", row.names = 1, quote = "")

ta_obj <- CreateSeuratObject(counts = expr_mtx, meta.data = metadata)

umap_coords <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/NG-2022/ascending_descending_human_aorta_umap.txt", skip = 1, header = TRUE, sep = "\t", row.names = 1)
umap_coords <- umap_coords[colnames(ta_obj), ]
colnames(umap_coords)[1:2] = c("UMAP_1", "UMAP_2")
ta_obj[['umap']] <- CreateDimReducObject(
  embeddings = as.matrix(umap_coords[, c("UMAP_1", "UMAP_2")]),  # Replace with actual column names
  key = "UMAP_",  # This is the prefix for UMAP coordinates
  assay = DefaultAssay(ta_obj)
)


cell_type_annotation <- c("03. VSMC II"="VSMC", "00. VSMC I"="VSMC", "04. Endothelial I"="EC", "01. Fibroblast I"="Fibroblast", "09. Endothelial II"="EC",
                         "05. Pericyte"="Pericyte", "13. ?"="Unknown", "02. Macrophage"="Macrophage", "08. Lymphatic Endothelial"="LEC", "06. Lymphocyte"="Lymphocyte",
                         "10. Neuron" = "Neuron", "12. ?"="Unknown", "07. Fibroblast II"="Fibroblast", "11. Mesothelial"="Mesothelial cell")
ta_obj$cell_type <- cell_type_annotation[as.character(ta_obj$cell_type_leiden)]
ta_obj$cell_type <- factor(ta_obj$cell_type, levels = c("EC", "LEC", "VSMC", "Pericyte", "Fibroblast", "Mesothelial cell", "Lymphocyte", "Macrophage", "Neuron", "Unknown"))

saveRDS(ta_obj, "/xdisk/mliang1/qqiu/project/others/test/data/NG-2022/ascending_descending_human_aorta.rds")







################################################################################
seurat_object = readRDS("/xdisk/mliang1/qqiu/project/others/test/cluster/huves1SN.cluster.rds")
ta_obj = readRDS("/xdisk/mliang1/qqiu/project/others/test/data/NG-2022/ascending_descending_human_aorta.rds")


### overall figures
cell_order <- levels(ta_obj$cell_type)
getPalette <- colorRampPalette(brewer.pal(12, "Paired"))
cell_col <- getPalette(length(cell_order))
names(cell_col) <- cell_order



DimPlot(seurat_object, group.by = "cell_type", label = T, repel = T) + labs(x = "UMAP 1", y = "UMAP 2", title = "Arterioles") + scale_color_manual(values = cell_col)
ggsave("/xdisk/mliang1/qqiu/project/others/test/figure/figure.1a.arterioles.umap.png", width=491/96, height=320/96, dpi=300)

DimPlot(ta_obj, group.by = "cell_type", label = T, repel = T) + labs(x = "UMAP 1", y = "UMAP 2", title = "Thoracic aorta") + scale_color_manual(values = cell_col)
ggsave("/xdisk/mliang1/qqiu/project/others/test/figure/figure.1a.thoracic_aorta.umap.png", width=491/96, height=320/96, dpi=300)



marker_list <- c("PECAM1", "EMCN", "PROX1", "ACTA2", "TAGLN", 
                "PDGFRB", "DCN", "PTPRC", "MRC1")
dot1 <- DotPlot(seurat_object, features=marker_list, group.by="cell_type") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="", title="Arterioles")
print(dot1)
ggsave("/xdisk/mliang1/qqiu/project/others/test/figure/figure.1b.arterioles.dotplot.png", width=507/96, height=349/96, dpi=300)

marker_list <- c("PECAM1", "EMCN", "PROX1", "ACTA2", "TAGLN", 
                "PDGFRB", "DCN", "WT1", "PTPRC", "MRC1", "MAPT")
dot2 <- DotPlot(ta_obj, features=marker_list, group.by="cell_type") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="", title = "Thoracic aorta")
print(dot2)
ggsave("/xdisk/mliang1/qqiu/project/others/test/figure/figure.1b.thoracic_aorta.dotplot.png", width=587/96, height=349/96, dpi=300)



### pericyte / VSMC to EC ratio
cell_count_arterioles <- table(seurat_object$cell_type)
cell_count_ta <- table(ta_obj$cell_type)

pe_ratio_arterioles <- cell_count_arterioles["Pericyte"] / cell_count_arterioles["EC"]
pe_ratio_ta <- cell_count_ta["Pericyte"] / cell_count_ta["EC"]
ve_ratio_arterioles <- cell_count_arterioles["VSMC"] / cell_count_arterioles["EC"]
ve_ratio_ta <- cell_count_ta["VSMC"] / cell_count_ta["EC"]

ratio_df <- data.frame(
  Dataset = c("Arterioles", "Thoracic Aorta",
              "Arterioles", "Thoracic Aorta"),
  Ratio = c(pe_ratio_arterioles, pe_ratio_ta,
            ve_ratio_arterioles, ve_ratio_ta),
  Comp = c("Pericyte/EC", "Pericyte/EC", 
           "VSMC/EC", "VSMC/EC")
)
ratio_df$Dataset <- factor(ratio_df$Dataset, levels = c("Thoracic Aorta", "Arterioles"))

ggplot(ratio_df, aes(x = Ratio, y = Dataset)) +
  geom_bar(stat = "identity", fill = "steelblue") + 
  theme_classic() + 
  labs(x = "Ratio",
       y = "") + 
  theme(axis.text.y = element_text(size = 12, color="black"),
        strip.text = element_text(size = 12, color="black")) +
  facet_wrap(~Comp, scales = "free_x")

ggsave("/xdisk/mliang1/qqiu/project/others/test/figure/figure.1c.ratio.barplot.png", width=364/96, height=271/96, dpi=300)



### pairwise correlation
pseudo_bulk_1 <- AverageExpression(seurat_object, return.seurat = F, group.by = c('cell_type'))$RNA
pseudo_bulk_2 <- AverageExpression(ta_obj, return.seurat = F, group.by = c('cell_type'))$RNA

common_genes <- intersect(rownames(pseudo_bulk_1), rownames(pseudo_bulk_2))
pseudo_bulk_1 <- pseudo_bulk_1[common_genes, ]
pseudo_bulk_2 <- pseudo_bulk_2[common_genes, ]

Idents(seurat_object) <- "cell_type"
degs_arterioles <- FindAllMarkers(seurat_object, group.by = "cell_type", assay = "RNA", only.pos = TRUE, min.pct = 0.25)
degs_arterioles <- rownames(degs_arterioles[degs_arterioles$p_val_adj < 0.05, ])
Idents(ta_obj) <- "cell_type"
degs_ta <- FindAllMarkers(ta_obj, group.by = "cell_type", only.pos = T)
degs_ta <- rownames(degs_ta[degs_ta$p_val_adj < 0.05, ])

common_deg <- intersect(common_genes, intersect(degs_arterioles, degs_ta))
write.table(common_deg, "/xdisk/mliang1/qqiu/project/others/test/cluster/cell_type_DEG.overlap.out", quote = F, row.names = F)


### pairwise correlation based on marker genes
common_deg <- read.table("/xdisk/mliang1/qqiu/project/others/test/cluster/cell_type_DEG.overlap.out", sep = "\t", header = T)
common_deg <- common_deg$x
pseudo_bulk_subset_1 <- pseudo_bulk_1[common_deg, ]
pseudo_bulk_subset_2 <- pseudo_bulk_2[common_deg, ]
all(rownames(pseudo_bulk_subset_1)==rownames(pseudo_bulk_subset_2))

cor_matrix <- cor(pseudo_bulk_subset_1, pseudo_bulk_subset_2, method = "spearman")

cor_matrix_melted <- melt(cor_matrix)
cor_matrix_melted$Var1 <- factor(cor_matrix_melted$Var1, levels = rev(levels(factor(cor_matrix_melted$Var1))))
ggplot(cor_matrix_melted, aes(Var2, Var1, fill = value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.3, 
                       limits = c(min(cor_matrix_melted$value), max(cor_matrix_melted$value))) +
  theme_minimal() +
  labs(title = "Expression Correlation", 
       x = "Thoracic Aorta", 
       y = "Arterioles", 
       fill = "Correlation") + 
  theme(axis.text = element_text(color = "black", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("/xdisk/mliang1/qqiu/project/others/test/figure/figure.1d.expression_cor-spearman.png", width=432/96, height=266/96, dpi=300)




