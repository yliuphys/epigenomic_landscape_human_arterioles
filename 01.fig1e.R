library(dplyr)
library(Seurat)
library(MAGMA.Celltyping)
library(SingleCellExperiment)

### .Renviron
# gh::gh_token()

storage_dir <- "/xdisk/mliang1/qqiu/project/others/test/MAGMA"



# munge GWAS
input_file <- c(

  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90104537_buildGRCh37.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90000065_buildGRCh37.tsv",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90000066_buildGRCh37.tsv",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90000063_buildGRCh37.tsv",
  
  "/xdisk/mliang1/qqiu/data/HT-GWAS/31217584-GCST008036-EFO_0000537-build37.f.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/31217584-GCST008044-EFO_0006335-build37.f.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/31217584-GCST008029-EFO_0006336-build37.f.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/29531354-GCST005841-EFO_1001504.h.tsv.gz",
  
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90310294.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90310295.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90310296.tsv.gz",
  
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90018832_buildGRCh37.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90435706.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90079903_buildGRCh38.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90435714.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90018890_buildGRCh37.tsv.gz"
  
)


for(i in input_file){
  
  gwas_sumstats_path <- i
  save_path <- gsub("vcf|txt.gz|tsv.gz|tsv|zip|txt", "formatted.tsv.gz", gwas_sumstats_path)

  tryCatch({
    path_formatted <- MungeSumstats::format_sumstats(path = gwas_sumstats_path,
                                                    save_path = save_path,
                                                    ref_genome = "GRCh38")
  }, error = function(e) {
    message(paste("Error processing file:", gwas_sumstats_path))
    message("Error message:", e$message)
  })

}


# map snps to genes
input_file = c(

  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90104537_buildGRCh37.formatted.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90000065_buildGRCh37.formatted.tsv",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90000066_buildGRCh37.formatted.tsv",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90000063_buildGRCh37.formatted.tsv",
  
  "/xdisk/mliang1/qqiu/data/HT-GWAS/31217584-GCST008036-EFO_0000537-build37.f.formatted.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/31217584-GCST008044-EFO_0006335-build37.f.formatted.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/31217584-GCST008029-EFO_0006336-build37.f.formatted.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/29531354-GCST005841-EFO_1001504.h.formatted.tsv.gz",
  
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90310294.formatted.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90310295.formatted.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90310296.formatted.tsv.gz",
  
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90018832_buildGRCh37.formatted.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90435706.formatted.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90079903_buildGRCh38.formatted.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90435714.formatted.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90018890_buildGRCh37.formatted.tsv.gz"
  
  )

input_df <- data.frame(input_file = input_file,
                   genome_ref = c("eur", "eur", "eur", "eur", "eur"),
                   N = c(585264, 388955, 387930, 398198, 660791))

for(i in 1:nrow(input_df)){
  
  gwas_sumstats_path <- input_df$input_file[i]
  genome_ref <- input_df$genome_ref[i]
  genome_ref_path <- paste0("/xdisk/mliang1/qqiu/reference/MAGMA.Celltyping/g1000_", genome_ref,"/g1000_", genome_ref)
  N <- input_df$N[i]
  ### need github token, set in .Renviron
  genesOutPath_intelligence <- MAGMA.Celltyping::map_snps_to_genes(
    path_formatted = gwas_sumstats_path,
    genome_ref_path = genome_ref_path,
    population = genome_ref,
    N = N,
    genome_build = "GRCh38")
  
}





# cell type dataset
input_file <- c(
  "/xdisk/mliang1/qqiu/project/others/test/cluster/huves1SN.cluster.rds",
  "/xdisk/mliang1/qqiu/project/others/test/data/NG-2022/ascending_descending_human_aorta.rds"
)

for(i in input_file){

  dataset <- gsub("\\..*ds", "", basename(i), perl = T)
  seurat_object <- readRDS(i)
  sce <- as.SingleCellExperiment(seurat_object, assay="RNA")
  sce <- fix_bad_mgi_symbols(sce)

  sce_dropped <- drop_uninformative_genes(exp=sce,
                                         input_species = "human",
                                         convert_orths = T,
                                         level2annot=sce$cell_type)

  annotLevels <- list(level1class=as.character(sce$cell_type)
                     )

  ctd <- generate_celltype_data(exp=sce_dropped,
                               annotLevels=annotLevels,
                               groupName=dataset,
                               savePath=storage_dir)

}








# https://neurogenomics.github.io/MAGMA_Celltyping/articles/full_workflow.html

input_file <- c(
  "/xdisk/mliang1/qqiu/project/others/test/MAGMA/ctd_huves1SN.rda",
  "/xdisk/mliang1/qqiu/project/others/test/MAGMA/ctd_ascending_descending_human_aorta.rda"
)


magma_dirs_list <- c(

  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90104537_buildGRCh37.formatted.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90000065_buildGRCh37.formatted.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90000066_buildGRCh37.formatted.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90000063_buildGRCh37.formatted.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/31217584-GCST008029-EFO_0006336-build37.f.formatted.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/31217584-GCST008044-EFO_0006335-build37.f.formatted.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/31217584-GCST008036-EFO_0000537-build37.f.formatted.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/29531354-GCST005841-EFO_1001504.h.formatted.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90310294.formatted.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90310295.formatted.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90310296.formatted.tsv.gz",
  
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90018832_buildGRCh37.formatted.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90435706.formatted.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90079903_buildGRCh38.formatted.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90435714.formatted.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90018890_buildGRCh37.formatted.tsv.gz"
  
)

for(magma_dirs in magma_dirs_list){

  result_merge <- c()
  for( i in input_file ){

    dataset <- gsub("ctd_(.+).rda", "\\1", basename(i), perl = T)
    ctd <- load_rdata(i)

    tryCatch({
      
      MAGMA_results <- MAGMA.Celltyping::celltype_associations_pipeline(
        magma_dirs = magma_dirs,
        ctd = ctd,
        force_new = T,
        ctd_name = dataset,
        run_linear = TRUE,
        run_top10 = TRUE,
        save_dir = storage_dir)
      
      merged_results <- MAGMA.Celltyping::merge_results(
        MAGMA_results = MAGMA_results)
      
      result_merge <- rbind(result_merge, merged_results)
      
    }, error = function(e) {
      message(paste("Error processing file:", i))
    })
    

  }

  outfile <- paste0("/xdisk/mliang1/qqiu/project/others/test/MAGMA/", gsub("formatted.tsv.gz", "", basename(magma_dirs)), "res.out")
  write.table(result_merge, outfile, row.names=F, col.names=T, sep='\t')

}



# other visualization method: https://star-protocols.cell.com/protocols/1392
base_font_size <- 12
theme_set(theme_classic(base_size = base_font_size))
setwd("/xdisk/mliang1/qqiu/project/others/test/MAGMA")

trait_use <- read.table("/xdisk/mliang1/qqiu/project/others/test/MAGMA/magma.trait_use.txt", header=T, sep='\t')

trait_use <- trait_use[! grepl("Wojcik|Haas", trait_use$Name), ]
result_file <- paste0("/xdisk/mliang1/qqiu/project/others/test/MAGMA/", trait_use$Files, ".res.out")

result_merge <- c()
for(i in result_file){

  result_tmp <- read.table(i, header = T, sep = '\t')
  
  file <- gsub(".res.out", "", basename(i))
  trait <- trait_use[trait_use$Files==file, ]$Name

  result_tmp$tissue <- ifelse(grepl("huves1SN", result_tmp$analysis_name), "Arterioles", "Thoracic aorta")
  result_tmp$tissue <- factor(result_tmp$tissue, levels=c("Arterioles", "Thoracic aorta"))
  result_tmp$trait <- factor(trait, levels=trait_use$Name)
  
  result_merge <- rbind(result_merge, result_tmp)

}

result_merge$study <- gsub(".*- ", "", result_merge$trait)
result_merge$trait <- factor(gsub(" -.*", "", result_merge$trait), 
levels <- c("Diastolic BP","Systolic BP","Pulse pressure","Hypertension","Stroke",
           "Diabetic nephropathy","Diabetic retinopathy","Peripheral Arterial Disease"))
result_merge[result_merge$Celltype_id=="PC", ]$Celltype_id <- "Pericyte"
result_merge[result_merge$Celltype_id=="FIB", ]$Celltype_id <- "Fibroblast"
p <- ggplot(result_merge[result_merge$EnrichmentMode=="Top 10%" & result_merge$FDR<0.05, ], # Linear
           aes(x = study, y = Celltype_id,
               size = BETA, fill = -1 * log10p)) +
  scale_fill_gradient(low = "white", high = "darkred") +
  scale_y_discrete(limits=rev) +
  geom_point(shape = 21) +
  theme(
    panel.grid.major.y = element_blank(), 
    panel.spacing.y=unit(0.2, "lines"),
    panel.spacing.x=unit(0.1, "lines"),
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black', size=12),
    axis.text.x = element_blank(),
    strip.text = element_text(colour = 'black', size=12)
  ) +
  labs(fill="-log10(p-value)", y="", x="") +
  facet_nested(tissue ~ trait + study,
             scales = "free", space = "free")
print(p)
ggsave("/xdisk/mliang1/qqiu/project/others/test/figure/fig1e.snRNA.MAGMA.png", width=754/96, height=390/96, dpi=300)




