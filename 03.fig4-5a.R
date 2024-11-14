library(dplyr)
library(reshape2)
library(GenomicRanges)
library(clusterProfiler)
library(org.Hs.eg.db)
library(patchwork)



################################################################################
#### preprocess
### load GWAS file and process

gwas_folder = "/xdisk/mliang1/qqiu/data/HT-GWAS/"
gwas_list = c("gwas_catalog_systolic_bp.tsv", "gwas_catalog_diastolic_bp.tsv",
              "gwas_catalog_pulse_pressure.tsv", "gwas_catalog_essential_hypertension.tsv",
              "gwas_catalog_stroke.tsv", "gwas_catalog_diabetic_nephropathy.tsv",
              "gwas_catalog_diabetic_retinopathy.tsv", "gwas_catalog_peripheral_arterial_disease.tsv")

gwas_folder = "/xdisk/mliang1/qqiu/data/HT-GWAS/control_traits/"
gwas_list = c("gwas_catalog_birth_weight.tsv", "gwas_catalog_bone_density.tsv",
              "gwas_catalog_breast_cancer.tsv", "gwas_catalog_colorectal_cancer.tsv",
              "gwas_catalog_alzheimer_disease.tsv", "gwas_catalog_parkinson_disease.tsv",
              "gwas_catalog_multiple_sclerosis.tsv", "gwas_catalog_rheumatoid_arthritis.tsv",
              "gwas_catalog_asthma.tsv", "gwas_catalog_psoriasis.tsv", "gwas_catalog_copd.tsv")

SNPS = c()
for(gwas_file in gwas_list){

  output = gsub("tsv", "processed.txt", gwas_file)
  trait = gsub("gwas_catalog_|\\.tsv", "", gwas_file)

  dat = read.table(paste0(gwas_folder, gwas_file), header=T, sep='\t', quote = "", fill = T)
  dat = dat[dat$P.VALUE<5e-8, ]

  snp_pos = data.frame(SNP=dat$SNPS, POS = paste0("chr", dat$CHR_ID, "_", dat$CHR_POS),
                       P.VALUE=dat$P.VALUE, BETA=dat$OR.or.BETA, TRAIT=trait, TYPE=dat$CONTEXT)

  snp_pos = snp_pos[(grepl("^rs", snp_pos$SNP)) & !(grepl("[,|;|x]", snp_pos$SNP) | grepl("NA", snp_pos$POS)), ]

  snp_pos_filtered <- snp_pos %>%
    group_by(SNP) %>%
    filter(P.VALUE == min(P.VALUE)) %>%
    ungroup()

  SNPS = rbind(SNPS, snp_pos_filtered)

}

write.table(SNPS, "/xdisk/mliang1/qqiu/project/others/test/data/gwas_catalog.processed.txt", col.names = T, row.names = F, quote = F, sep = '\t')




gwas = read.table("/xdisk/mliang1/qqiu/project/others/test/data/gwas_catalog.processed.txt", header = T, sep = "\t")
gwas <- gwas %>%
  mutate(CHROM = sub("_.*", "", POS),
         POS = as.numeric(sub(".*_", "", POS)))  %>%
  filter(!is.na(POS))
gwas_ranges <- makeGRangesFromDataFrame(gwas,
                                        seqnames.field = "CHROM",
                                        start.field = "POS",
                                        end.field = "POS",  # SNP is a single-point
                                        keep.extra.columns = TRUE)

# mcols(gwas_ranges)$SNP_CHR <- as.character(seqnames(gwas_ranges))
# mcols(gwas_ranges)$SNP_POS <- start(gwas_ranges)

regulatory_ref = read.table("/xdisk/mliang1/qqiu/project/others/test/data/regulatory.reference.out", header = F, sep = "\t")
regulatory_ref = regulatory_ref[regulatory_ref$V3!="promoter", ]
regulatory_ref = regulatory_ref[, c(1, 3:5)]
regulatory_ref$V1 = paste0("chr", regulatory_ref$V1)

regulatory_ranges <- makeGRangesFromDataFrame(
  regulatory_ref,
  seqnames.field = "V1",
  start.field = "V4",
  end.field = "V5",
  keep.extra.columns = TRUE
)


canonical_chromosomes <- c(paste0("chr", 1:22), "chrX", "chrY")
ref <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/Homo_sapiens.GRCh38.112.transcript.gtf", header = F, sep = '\t')

transcript_ranges <- GRanges(
  seqnames = Rle(paste0("chr", ref$V1)),
  ranges = IRanges(start = ref$V4, end = ref$V5),
  strand = Rle(ref$V7)
)
promoter_ranges <- promoters(transcript_ranges, upstream = 1000, downstream = 500)
mcols(promoter_ranges)$V3 = "promoter"
promoter_ranges <- promoter_ranges[seqnames(promoter_ranges) %in% canonical_chromosomes]
promoter_ranges <- keepSeqlevels(promoter_ranges, canonical_chromosomes, pruning.mode = "coarse")

regulatory_ranges <- c(regulatory_ranges, promoter_ranges)

colnames(mcols(regulatory_ranges)) = "REG_ANNO"
# mcols(regulatory_ranges)$REG_CHR <- as.character(seqnames(regulatory_ranges))
# mcols(regulatory_ranges)$REG_START <- start(regulatory_ranges)
# mcols(regulatory_ranges)$REG_END <- end(regulatory_ranges)

overlaps <- findOverlaps(gwas_ranges, regulatory_ranges)

gwas_reg <- cbind(as.data.frame(mcols(gwas_ranges))[queryHits(overlaps), ],
                  data.frame(REG_ANNO=mcols(regulatory_ranges)[subjectHits(overlaps), ]))

merge_list = intersect(colnames(gwas), colnames(gwas_reg))
gwas_merge = merge(gwas, gwas_reg, by = merge_list, all.x=T)

write.table(gwas_merge, "/xdisk/mliang1/qqiu/project/others/test/data/gwas_catalog.processed.txt", col.names = T, row.names = F, quote = F, sep = '\t')


### load expr file and process
library(DESeq2)

count_data <- read.csv("/xdisk/mliang1/qqiu/project/others/test/data/vessel.gene_count_matrix.cln.csv", header = TRUE, row.names = 1, sep = ",")
rownames(count_data) <- gsub("\\..*", "", sapply(rownames(count_data), function(x) strsplit(x, "\\|")[[1]][1]))
condition <- factor(rep(c("EC", "NA", "VR"), 4))

dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = data.frame(condition),
                              design = ~ 1)

dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)
head(normalized_counts)

write.csv(normalized_counts, file = "/xdisk/mliang1/qqiu/project/others/test/data/vessel.gene.normalized_counts.csv")



### load methylation file and process
canonical = read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene_tran.canonical.GRCH38.out", header = T, sep = "\t")
canonical = canonical[!is.na(canonical$Ensembl.Canonical),]

meth = read.table("/xdisk/mliang1/qqiu/project/others/test/data/vessel.RRBS.dat.cln.promoter.methyl.cln", header = T)
meth$transcript = gsub("\\..*", "", meth$transcript)
meth$gene = gsub("\\..*", "", meth$gene)
meth = meth[paste(meth$transcript, meth$gene) %in%
              paste(canonical$Transcript.stable.ID, canonical$Gene.stable.ID), ]
rownames(meth) = meth$gene
meth = meth[, 10:21]
write.csv(meth, file = "/xdisk/mliang1/qqiu/project/others/test/data/vessel.RRBS.dat.cln.promoter.methyl.out")











################################################################################
#### trait level enrichment

### functions
standardize_column_names <- function(df) {
  df %>%
    rename_with(~"BIN1_CHR",  matches("bait_chr|bin1_chr|chr1|BIN1_CHR")) %>%
    rename_with(~"BIN1_START", matches("bait_start|bin1_start|start1|BIN1_START")) %>%
    rename_with(~"BIN1_END", matches("bait_end|bin1_end|end1|BIN1_END")) %>%
    
    rename_with(~"BIN2_CHR",  matches("otherEnd_chr|bin2_chr|chr2|BIN2_CHR")) %>%
    rename_with(~"BIN2_START", matches("otherEnd_start|bin2_start|start2|BIN2_START")) %>%
    rename_with(~"BIN2_END", matches("otherEnd_end|bin2_end|end2|BIN2_END"))
}

calculate_p_value <- function(covered_snps, total_snps, control_mean) {
  test_result <- binom.test(covered_snps, total_snps, control_mean, "greater")
  
  return(test_result$p.value)
}



### load all files
mc_NV <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/micro-c.merged_NV_8kb_loops.tsv", header = T)
mc_VR <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/micro-c.merged_VR_8kb_loops.tsv", header = T)

cm_NV <- read.csv("/xdisk/mliang1/qqiu/project/others/test/data/NVCMC_newpara_20kb_intra.filtered.csv", header = T)
cm_VR <- read.csv("/xdisk/mliang1/qqiu/project/others/test/data/VRCMC_newpara_20kb_intra.filtered.csv", header = T)

gwas <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/gwas_catalog.processed.txt", header = T, sep = "\t")



### trait level enrichment
trait_list <- c("systolic_bp", "diastolic_bp", "pulse_pressure", "essential_hypertension", "stroke", 
               "diabetic_nephropathy", "diabetic_retinopathy", "peripheral_arterial_disease")

input_df <- data.frame(input = c("mc_NV", "mc_VR", "cm_NV", "cm_VR"),
                      title = c("Micro-C (arteriole)", "Micro-C (EDA)", "Capture Micro-C (arteriole)", "Capture Micro-C (EDA)"))


for( i in 1:nrow(input_df)){
  
  input <- get(input_df[i,]$input)
  outfile <- paste0(input_df[i,]$input, ".traits.bar_plot.png")
  title <- input_df[i,]$title
  fill_col <- ifelse(grepl("arteriole", title), "skyblue2", "salmon2")
  
  input <- standardize_column_names(input)
  
  mc_interactions <- makeGRangesFromDataFrame(
    input,
    seqnames.field = "BIN1_CHR",
    start.field = "BIN1_START",
    end.field = "BIN1_END"
  )
  
  mc_interactions_bin2 <- makeGRangesFromDataFrame(
    input,
    seqnames.field = "BIN2_CHR",
    start.field = "BIN2_START",
    end.field = "BIN2_END"
  )
  
  traits <- unique(gwas_data$TRAIT)
  results <- list()
  for (trait in traits) {
    trait_gwas <- gwas_data %>% filter(TRAIT == trait)
    gwas_snps <- makeGRangesFromDataFrame(
      trait_gwas,
      seqnames.field = "CHROM",
      start.field = "POS_START",
      end.field = "POS_START",
      keep.extra.columns = TRUE
    )
    
    overlaps_bin1 <- findOverlaps(gwas_snps, mc_interactions)
    overlaps_bin2 <- findOverlaps(gwas_snps, mc_interactions_bin2)
    
    covered_snps <- length(unique(c(queryHits(overlaps_bin1), queryHits(overlaps_bin2))))
    total_snps <- length(unique(gwas_snps$SNP))
    coverage_rate <- covered_snps / total_snps
    
    results[[trait]] <- data.frame(
      Trait = trait,
      Covered_SNPs = covered_snps,
      Total_SNPs = total_snps,
      Coverage_Rate = coverage_rate
    )
  }
  coverage_results <- bind_rows(results)
  
  # print(coverage_results)
  
  coverage_results %>%
    filter(Coverage_Rate > 0) %>%
    mutate(Group = ifelse(Trait %in% trait_list, "Interested", "Control")) -> processed_data
  
  mean_control_coverage <- processed_data %>%
    filter(Group == "Control") %>%
    summarize(mean_coverage = mean(Coverage_Rate)) %>%
    pull(mean_coverage)
  
  processed_data <- processed_data %>% filter(Group == "Interested") %>% 
    rowwise() %>%
    mutate(
      p_value = calculate_p_value(Covered_SNPs, Total_SNPs, mean_control_coverage),
      significance = case_when(
        p_value < 0.001 ~ "***",  # Highly significant
        p_value < 0.01  ~ "**",   # Significant
        p_value < 0.05  ~ "*",    # Marginally significant
        TRUE            ~ ""      # Not significant
      )
    ) %>%
    ungroup()
  
  p <- ggplot(processed_data,
         aes(x = fct_reorder(Trait, Coverage_Rate), y = Coverage_Rate)) +
    geom_bar(stat = "identity", fill=fill_col) +
    geom_hline(yintercept = mean_control_coverage, linetype = "dashed", color = "red", size = 1) +
    geom_text(aes(label = significance, y = Coverage_Rate+0.1*max(Coverage_Rate)), 
              color = "black", size = 5, hjust = 1) +
    theme_classic() + 
    theme(axis.text.y = element_text(color="black", size=12)) +
    labs(title = title,
         y = "Percent of SNPs covered by interactions", x = "") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
    coord_flip()
  
  print(p)
  
  ggsave(paste0("/xdisk/mliang1/qqiu/project/others/test/figure/fig4a.", outfile), width=528/96, height=227/96, dpi=300)
  
}



### snp per 1m across datasets
data <- data.frame(
  Method = c("Tibial artery (HiC, 5kb)", "Aorta 1 (HiC, 2kb)", "Aorta 2 (HiC, 2kb)", 
             "Arteriole (Micro-C, 4kb)", "Arteriole (Micro-C, 8kb)", 
             "EDA (Micro-C, 4kb)", "EDA (Micro-C, 8kb)"),
  SNP_per_1M_contact = c(0.917654551, 1.010804321, 1.066421048, 
                         1.58556, 1.69930, 1.33599, 1.19008)
)
data$Method <- factor(data$Method, levels = data$Method)
hic_baseline <- mean(data$SNP_per_1M_contact[grep("HiC", data$Method)]) 

p <- ggplot(data, aes(x = Method, y = SNP_per_1M_contact)) +
  geom_bar(stat = "identity", fill = "grey") +
  geom_hline(yintercept = hic_baseline, linetype = "dashed", color = "red", size = 1) +
  annotate("text", x = 1, y = 1.4, label = "Large arteries", color = "red", size = 4) +
  theme_classic() +
  theme(axis.text = element_text(color="black", size=10),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "Number of BP SNPs per 1Mbp\nloop contact region") +
  coord_flip()
print(p)
ggsave("/xdisk/mliang1/qqiu/project/others/test/figure/fig4b.snp_per_1m.png", width=410/96, height=219/96, dpi=300)






################################################################################
#### function annotation
### functions
standardize_column_names <- function(df) {
  df %>%
    rename_with(~"BIN1_CHR",  matches("bait_chr|bin1_chr|chr1|BIN1_CHR")) %>%
    rename_with(~"BIN1_START", matches("bait_start|bin1_start|start1|BIN1_START")) %>%
    rename_with(~"BIN1_END", matches("bait_end|bin1_end|end1|BIN1_END")) %>%
    
    rename_with(~"BIN2_CHR",  matches("otherEnd_chr|bin2_chr|chr2|BIN2_CHR")) %>%
    rename_with(~"BIN2_START", matches("otherEnd_start|bin2_start|start2|BIN2_START")) %>%
    rename_with(~"BIN2_END", matches("otherEnd_end|bin2_end|end2|BIN2_END"))
}

### load reference files
canonical_chromosomes <- c(1:22, "X", "Y")

ref <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/Homo_sapiens.GRCh38.112.transcript.gtf", sep = '\t', header=F)
ref <- ref[ref$V1 %in% canonical_chromosomes, ]
ref$transcript_id <- str_extract(ref$V9, "transcript_id [^;]+") %>%
  gsub("transcript_id ", "", .)
ref$gene_id <- str_extract(ref$V9, "gene_id [^;]+") %>%
  gsub("gene_id ", "", .)
ref$gene_name <- str_extract(ref$V9, "gene_name [^;]+") %>%
  gsub("gene_name ", "", .)
ref$gene_type <- str_extract(ref$V9, "gene_biotype [^;]+") %>%
  gsub("gene_biotype ", "", .)
# rownames(ref) = ref$Gene.stable.ID

gene_ranges <- GRanges(
  seqnames = Rle(paste0("chr", ref$V1)),
  ranges = IRanges(start = ref$V4, end = ref$V5),
  strand = Rle(ref$V7),
  transcript_id = ref$transcript_id,
  gene_id = ref$gene_id,
  gene_name = ref$gene_name,
  gene_type = ref$gene_type,
  gene_chr = paste0("chr", ref$V1),
  gene_start = ref$V4,
  gene_end = ref$V5,
  gene_strand = ref$V7,
  overlap_type = ""
)

promoter_ranges <- promoters(gene_ranges, upstream = 1000, downstream = 500)
mcols(promoter_ranges)$overlap_type = "promoter"

gene_ranges <- c(gene_ranges, promoter_ranges)

### load gwas
trait_list <- c("systolic_bp", "diastolic_bp", "pulse_pressure", "essential_hypertension", "stroke", 
               "diabetic_nephropathy", "diabetic_retinopathy", "peripheral_arterial_disease")
gwas <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/gwas_catalog.processed.txt", header = T, sep = "\t")
gwas <- gwas[gwas$TRAIT %in% trait_list, ]

gwas_ranges <- makeGRangesFromDataFrame(gwas, 
                                        seqnames.field = "CHROM", 
                                        start.field = "POS", 
                                        end.field = "POS",  # SNP is a single-point
                                        keep.extra.columns = TRUE)
mcols(gwas_ranges)$SNP_CHR <- as.character(seqnames(gwas_ranges)) 
mcols(gwas_ranges)$SNP_POS <- start(gwas_ranges)    

### load other files
mc_NV <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/micro-c.merged_NV_8kb_loops.tsv", header = T)
mc_VR <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/micro-c.merged_VR_8kb_loops.tsv", header = T)

cm_NV <- read.csv("/xdisk/mliang1/qqiu/project/others/test/data/NVCMC_newpara_20kb_intra.filtered.csv", header = T)
cm_VR <- read.csv("/xdisk/mliang1/qqiu/project/others/test/data/VRCMC_newpara_20kb_intra.filtered.csv", header = T)

cm_NV <- cm_NV[, setdiff(colnames(cm_NV), c("start", "end"))]
cm_VR <- cm_VR[, setdiff(colnames(cm_VR), c("start", "end"))]

TF_list <- read.table("/xdisk/mliang1/qqiu/reference/Homo_sapiens_TF.txt", header = T, sep = "\t")

### generate SNP-gene interaction table
input_list <- c("mc_NV", "mc_VR", "cm_NV", "cm_VR")

for( i in input_list){
  
  input <- get(i)
  outfile <- paste0(i, ".snp_gene.out")
  
  input <- standardize_column_names(input) %>%
    filter(BIN1_CHR==BIN2_CHR) %>%
    mutate(mc_id = paste0(i, "_", row_number()))
  
  mc_ranges_bin1 <- makeGRangesFromDataFrame(input, 
                                             seqnames.field = "BIN1_CHR", 
                                             start.field = "BIN1_START", 
                                             end.field = "BIN1_END", 
                                             keep.extra.columns = TRUE)
  mcols(mc_ranges_bin1)$BIN_ID = "BIN1"
  
  mc_ranges_bin2 <- makeGRangesFromDataFrame(input, 
                                             seqnames.field = "BIN2_CHR", 
                                             start.field = "BIN2_START", 
                                             end.field = "BIN2_END", 
                                             keep.extra.columns = TRUE)
  mcols(mc_ranges_bin2)$BIN_ID = "BIN2"
  
  mc_ranges <- c(mc_ranges_bin1, mc_ranges_bin2)
  mcols(mc_ranges)$BIN1_CHR <- ifelse(is.na(mcols(mc_ranges)$BIN1_CHR), as.character(seqnames(mc_ranges)), mcols(mc_ranges)$BIN1_CHR)
  mcols(mc_ranges)$BIN1_START <- ifelse(is.na(mcols(mc_ranges)$BIN1_START), start(mc_ranges), mcols(mc_ranges)$BIN1_START)
  mcols(mc_ranges)$BIN1_END <- ifelse(is.na(mcols(mc_ranges)$BIN1_END), end(mc_ranges), mcols(mc_ranges)$BIN1_END)
  mcols(mc_ranges)$BIN2_CHR <- ifelse(is.na(mcols(mc_ranges)$BIN2_CHR), as.character(seqnames(mc_ranges)), mcols(mc_ranges)$BIN2_CHR)
  mcols(mc_ranges)$BIN2_START <- ifelse(is.na(mcols(mc_ranges)$BIN2_START), start(mc_ranges), mcols(mc_ranges)$BIN2_START)
  mcols(mc_ranges)$BIN2_END <- ifelse(is.na(mcols(mc_ranges)$BIN2_END), end(mc_ranges), mcols(mc_ranges)$BIN2_END)
  
  gwas_mc_overlaps <- findOverlaps(mc_ranges, gwas_ranges)
  mc_hits_gwas <- mc_ranges[queryHits(gwas_mc_overlaps)]
  snp_hits <- gwas_ranges[subjectHits(gwas_mc_overlaps)]
  mc_gwas <- cbind(mcols(mc_hits_gwas), mcols(snp_hits))
  colnames(mc_gwas)[colnames(mc_gwas)=="BIN_ID"] = "BIN_ID_SNP"
  
  gene_mc_overlaps <- findOverlaps(mc_ranges, gene_ranges)
  mc_hits_gene <- mc_ranges[queryHits(gene_mc_overlaps)]
  gene_hits <- gene_ranges[subjectHits(gene_mc_overlaps )]
  mc_gene <- cbind(mcols(mc_hits_gene), mcols(gene_hits))
  colnames(mc_gene)[colnames(mc_gene)=="BIN_ID"] = "BIN_ID_GENE"
  
  by_list = intersect(colnames(mc_gwas), colnames(mc_gene))
  snp_gene = merge(mc_gwas, mc_gene, by=by_list)
  snp_gene = snp_gene[snp_gene$BIN_ID_SNP!=snp_gene$BIN_ID_GENE, ]
  snp_gene = snp_gene %>% as.data.frame() %>%
    filter(SNP_CHR==gene_chr) %>%
    mutate(TSS = ifelse(gene_strand == 1, gene_start, gene_end),
           distance = abs(SNP_POS - TSS)) %>%
    group_by(across(-c(transcript_id, gene_chr, gene_start, gene_end, gene_strand, TSS, distance, overlap_type))) %>%
    mutate(has_promoter_overlap = any(overlap_type == "promoter")) %>%
    filter(ifelse(has_promoter_overlap, overlap_type == "promoter", TRUE)) %>%
    slice_min(order_by = distance, with_ties = F) %>%
    ungroup() %>%
    mutate(distance_category = cut(distance,
                                   breaks = c(0, 10000, 50000, 100000, 200000, 500000, Inf),
                                   labels = c("<10 kb", "10-50 kb", "50-100 kb", "100-200 kb", "200-500 kb", "> 500kb"),
                                   include.lowest = TRUE)) %>%
    mutate(is_TF = ifelse(gene_id %in% TF_list$Ensembl, TRUE, FALSE)) %>%
    unique()
  
  write.table(snp_gene, paste0("/xdisk/mliang1/qqiu/project/others/test/integration/", outfile), col.names = T, row.names = F, sep = "\t", quote = F)

}



################################################################################
### overlap SNP-gene pair with those in GWAS Catalog and GTEx-eQTL
#### preprocess GWAS snp-gene pairs
ensembl <- read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.GRCH38.out", sep = '\t', header=T)

gwas_folder <- "/xdisk/mliang1/qqiu/data/HT-GWAS/"
gwas_list <- c("gwas_catalog_systolic_bp.tsv", "gwas_catalog_diastolic_bp.tsv",
              "gwas_catalog_pulse_pressure.tsv", "gwas_catalog_essential_hypertension.tsv",
              "gwas_catalog_stroke.tsv", "gwas_catalog_diabetic_nephropathy.tsv",
              "gwas_catalog_diabetic_retinopathy.tsv", "gwas_catalog_peripheral_arterial_disease.tsv")
snp_gene_pairs_all <- c()
for(gwas_file in gwas_list){
  
  output <- gsub("tsv", "processed.txt", gwas_file)
  trait <- gsub("gwas_catalog_|\\.tsv", "", gwas_file)
  
  dat <- read.table(paste0(gwas_folder, gwas_file), header=T, sep='\t', quote = "", fill = T)
  dat_filter = dat %>% filter(P.VALUE<5e-8) %>%
    filter(grepl("^rs", SNPS), !(grepl("[,|;|x]", SNPS) | grepl("NA", CHR_POS) | grepl("NA", CHR_ID))) %>%
    dplyr::select(SNPS, REPORTED.GENE.S., MAPPED_GENE, UPSTREAM_GENE_ID, DOWNSTREAM_GENE_ID, SNP_GENE_IDS)
  
  long_gene_data <- dat_filter %>%
    gather(key = "Gene_Type", value = "Gene", -SNPS) %>%
    filter(!is.na(Gene) & Gene != "") %>%
    separate_rows(Gene, sep = ",\\s*|\\s+-\\s+")
    
  ensembl_genes <- long_gene_data %>%
    filter(grepl("^ENSG", Gene))
  
  symbol_genes <- long_gene_data %>%
    filter(!grepl("^ENSG", Gene))
  
  symbol_to_ensembl <- symbol_genes %>%
    left_join(ensembl, by = c("Gene" = "Gene.name")) %>%
    filter(!is.na(Gene.stable.ID)) %>% 
    mutate(Gene_Name = Gene,
           Gene = Gene.stable.ID) %>%
    dplyr::select(SNPS, Gene_Type, Gene, Gene_Name)
  
  ensembl_to_symbol <- ensembl_genes %>%
    left_join(ensembl, by = c("Gene" = "Gene.stable.ID")) %>%
    mutate(Gene_Name = Gene.name) %>%  
    dplyr::select(SNPS, Gene_Type, Gene, Gene_Name)
  
  all_ensembl_genes <- bind_rows(symbol_to_ensembl, ensembl_to_symbol)
  
  snp_gene_pairs <- all_ensembl_genes %>%
    dplyr::select(SNPS, Gene_Type, Gene, Gene_Name) %>%
    distinct()
  
  snp_gene_pairs_all <- bind_rows(snp_gene_pairs_all, snp_gene_pairs)
  
}

snp_gene_pairs_aggregated <- snp_gene_pairs_all %>%
  distinct() %>%
  group_by(SNPS, Gene, Gene_Name) %>%
  summarize(Gene_Type = paste(unique(Gene_Type), collapse = ", "), .groups = 'drop') # Aggregating Gene_Type

write.table(snp_gene_pairs_aggregated, "/xdisk/mliang1/qqiu/project/others/test/data/gwas_catalog.arteriole_related.snp_gene.txt", col.names = T, row.names = F, quote = F, sep = '\t')



#### preprocess eQTL snp-gene pairs
process_tissue_data <- function(df, Tissue_Type) {
  df %>%
    mutate(Tissue_Type = Tissue_Type) %>%  
    dplyr::select(rs_id, gene_id, gene_name, Tissue_Type) %>%  
    mutate(gene_id = gsub("\\..*", "", gene_id)) %>%  
    rename(SNPS = rs_id, Gene = gene_id, Gene_Name = gene_name) %>%  
    distinct() 
}

trait_list <- c("systolic_bp", "diastolic_bp", "pulse_pressure", "essential_hypertension", "stroke", 
               "diabetic_nephropathy", "diabetic_retinopathy", "peripheral_arterial_disease")

gwas_snps <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/gwas_catalog.processed.txt", header = T, sep = "\t")
gwas_snps <- gwas_snps[gwas_snps$TRAIT %in% trait_list, ]

aorta <- read.table("/xdisk/mliang1/qqiu/data/GTEx/GTEx_Analysis_v8_eQTL/Artery_Aorta.v8.signif_variant_gene_pairs.processed.txt", sep = "\t", header = T)
tibial <- read.table("/xdisk/mliang1/qqiu/data/GTEx/GTEx_Analysis_v8_eQTL/Artery_Tibial.v8.signif_variant_gene_pairs.processed.txt", sep = "\t", header = T)
coronary <- read.table("/xdisk/mliang1/qqiu/data/GTEx/GTEx_Analysis_v8_eQTL/Artery_Coronary.v8.signif_variant_gene_pairs.processed.txt", sep = "\t", header = T)

aorta_processed <- process_tissue_data(aorta, "Artery_Aorta")
tibial_processed <- process_tissue_data(tibial, "Artery_Tibial")
coronary_processed <- process_tissue_data(coronary, "Artery_Coronary")

snp_gene_pairs_all <- bind_rows(aorta_processed, tibial_processed, coronary_processed)

snp_gene_pairs_aggregated <- snp_gene_pairs_all %>%
  filter(SNPS %in% gwas_snps$SNP) %>%
  distinct() %>%
  group_by(SNPS, Gene, Gene_Name) %>%
  summarize(Tissue_Type = paste(unique(Tissue_Type), collapse = ", "), .groups = 'drop')

write.table(snp_gene_pairs_aggregated, "/xdisk/mliang1/qqiu/project/others/test/data/eQTL.arteriole_related.snp_gene.txt", col.names = T, row.names = F, quote = F, sep = '\t')




#### preprocess chromatin-interaction based snp-gene pairs
process_inter_data <- function(df, Interaction_Type) {
  df %>%
    mutate(Interaction_Type = Interaction_Type) %>%  
    dplyr::select(SNP, gene_id, gene_name, Interaction_Type) %>%  
    rename(SNPS = SNP, Gene = gene_id, Gene_Name = gene_name) %>%  
    distinct() 
}

mc_NV <- read.table("/xdisk/mliang1/qqiu/project/others/test/integration/mc_NV.snp_gene.out", sep = "\t", header = T)
mc_VR <- read.table("/xdisk/mliang1/qqiu/project/others/test/integration/mc_VR.snp_gene.out", sep = "\t", header = T)
cm_NV <- read.table("/xdisk/mliang1/qqiu/project/others/test/integration/cm_NV.snp_gene.out", sep = "\t", header = T)
cm_VR <- read.table("/xdisk/mliang1/qqiu/project/others/test/integration/cm_VR.snp_gene.out", sep = "\t", header = T)

mc_NV_processed <- process_inter_data(mc_NV, "Micro_C_arteriole_8kb")
mc_VR_processed <- process_inter_data(mc_VR, "Micro_C_EDA_8kb")
cm_NV_processed <- process_inter_data(cm_NV, "Capture_Micro_C_arteriole_20kb")
cm_VR_processed <- process_inter_data(cm_VR, "Capture_Micro_C_EDA_20kb")

snp_gene_pairs_all <- combined_data <- bind_rows(mc_NV_processed, mc_VR_processed, cm_NV_processed, cm_VR_processed)

snp_gene_pairs_aggregated <- snp_gene_pairs_all %>%
  distinct() %>%
  group_by(SNPS, Gene, Gene_Name) %>%
  summarize(Interaction_Type = paste(unique(Interaction_Type), collapse = ", "), .groups = 'drop')

write.table(snp_gene_pairs_aggregated, "/xdisk/mliang1/qqiu/project/others/test/data/mc_cm.arteriole_related.snp_gene.txt", col.names = T, row.names = F, quote = F, sep = '\t')




### statistical overlap between eqtl and interaction snp-gene pairs
#### generate backround snp-gene pair (2m)
trait_list <- c("systolic_bp", "diastolic_bp", "pulse_pressure", "essential_hypertension", "stroke", 
                "diabetic_nephropathy", "diabetic_retinopathy", "peripheral_arterial_disease")

gwas_snps <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/gwas_catalog.processed.txt", header = T, sep = "\t")
gwas_snps <- gwas_snps[gwas_snps$TRAIT %in% trait_list, ]

ensembl <- read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.GRCH38.out", sep = '\t', header=T)


gwas_snps <- gwas_snps %>%
  mutate(POS = as.numeric(POS))

ensembl <- ensembl %>%
  mutate(Chromosome = paste0("chr", Chromosome.scaffold.name)) %>%
  rename(Gene_ID = Gene.stable.ID, Gene_start = Gene.start..bp., Gene_end = Gene.end..bp.) %>%
  mutate(Gene_start = as.numeric(Gene_start), Gene_end = as.numeric(Gene_end))

snp_gene_pairs <- gwas_snps %>%
  inner_join(ensembl, by = c("CHROM" = "Chromosome")) %>%
  filter(abs(POS - Gene_start) <= 2000000 | abs(POS - Gene_end) <= 2000000)  # SNP and gene within 2Mb

snp_gene_pairs <- snp_gene_pairs %>%
  mutate(distance_to_gene = pmin(abs(POS - Gene_start), abs(POS - Gene_end))) %>%
  dplyr::select(SNP, Gene_ID, Gene.name, CHROM, POS, Gene_start, Gene_end, distance_to_gene)

write.csv(snp_gene_pairs, "/xdisk/mliang1/qqiu/project/others/test/data/snp_gene_pairs_within_2mb.csv", row.names = FALSE)





snp_gene_pairs_bg <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/snp_gene_pairs_within_2mb.csv", header = T, sep = ",")
snp_gene_pairs_eqtl <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/eQTL.arteriole_related.snp_gene.txt", header = T, sep = "\t")
snp_gene_pairs_ci <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/mc_cm.arteriole_related.snp_gene.txt", header = T, sep = "\t")
snp_gene_pairs_gwas <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/gwas_catalog.arteriole_related.snp_gene.txt", header = T, sep = "\t")


snp_gene_pairs_bg <- snp_gene_pairs_bg %>%
  mutate(pair_id = paste(SNP, Gene_ID, sep = "_"),
         pair_symbol = paste(SNP, Gene.name, sep = "_"))

snp_gene_pairs_eqtl <- snp_gene_pairs_eqtl %>%
  mutate(pair_id = paste(SNPS, Gene, sep = "_"),
         pair_symbol = paste(SNPS, Gene_Name, sep = "_"))

snp_gene_pairs_ci <- snp_gene_pairs_ci %>%
  mutate(pair_id = paste(SNPS, Gene, sep = "_"),
         pair_symbol = paste(SNPS, Gene_Name, sep = "_"))

snp_gene_pairs_gwas <- snp_gene_pairs_gwas %>%
  mutate(pair_id = paste(SNPS, Gene, sep = "_"),
         pair_symbol = paste(SNPS, Gene_Name, sep = "_"))


length(unique(snp_gene_pairs_ci$pair_id))
sum(!(snp_gene_pairs_ci$pair_id %in% c(snp_gene_pairs_gwas$pair_id, snp_gene_pairs_eqtl$pair_id)))
# 10217/1684

length(unique(intersect(snp_gene_pairs_eqtl$pair_id, snp_gene_pairs_ci$pair_id)))
sum(intersect(snp_gene_pairs_eqtl$pair_id, snp_gene_pairs_ci$pair_id) %in% snp_gene_pairs_gwas$pair_id)
# 938/416

length(unique(intersect(snp_gene_pairs_eqtl$pair_symbol, snp_gene_pairs_ci$pair_symbol)))
sum(intersect(snp_gene_pairs_eqtl$pair_symbol, snp_gene_pairs_ci$pair_symbol) %in% snp_gene_pairs_gwas$pair_symbol)
# 899/416

snp_gene_pairs_bg_list <- unique(snp_gene_pairs_bg$pair_id)
snp_gene_pairs_eqtl_list <- unique(snp_gene_pairs_eqtl$pair_id)
snp_gene_pairs_ci_list <- unique(snp_gene_pairs_ci$pair_id)
overlap <- intersect(snp_gene_pairs_eqtl_list, snp_gene_pairs_ci_list)

eqtl_list_length <- length(snp_gene_pairs_eqtl_list)
ci_list_length <- length(snp_gene_pairs_ci_list)
bg_list_length <- length(snp_gene_pairs_bg_list)
overlap_length <- length(overlap)

p_value <- phyper(overlap_length - 1, ci_list_length, bg_list_length - ci_list_length, eqtl_list_length, lower.tail = FALSE)



only_eqtl <- eqtl_list_length - overlap_length
only_ci <- ci_list_length - overlap_length
neither <- bg_list_length - eqtl_list_length - ci_list_length + overlap_length
fisher_matrix <- matrix(c(overlap_length, only_eqtl, only_ci, neither), nrow = 2, byrow = TRUE)

fisher_result <- fisher.test(fisher_matrix)




library(ggVennDiagram)
venn_data <- list(snp_gene_pairs_eqtl_list, snp_gene_pairs_ci_list)
names(venn_data) <- c("eQTL", "Chromatin interaction")

ggVennDiagram(venn_data, label="both", label_alpha = 0) +
  scale_fill_gradientn(colors = c("#b3cde3", "#66c2a5", "#8da0cb")) +
  theme(legend.position = "None",) + 
  # labs(title = "Overlap of eQTL and chromatin interaction SNP-gene pairs\n(odds ratio = 11, p-value < 2e-16)") +
  coord_flip()
ggsave("/xdisk/mliang1/qqiu/project/others/test/figure/fig1c.eqtl_interaction.overlap.venn.png", width=349/96, height=224/96, dpi=300)




### additional SNP-associated genes and functional annotation
mc_NV <- read.table("/xdisk/mliang1/qqiu/project/others/test/integration/mc_NV.snp_gene.out", sep = "\t", header = T)
mc_VR <- read.table("/xdisk/mliang1/qqiu/project/others/test/integration/mc_VR.snp_gene.out", sep = "\t", header = T)
cm_NV <- read.table("/xdisk/mliang1/qqiu/project/others/test/integration/cm_NV.snp_gene.out", sep = "\t", header = T)
cm_VR <- read.table("/xdisk/mliang1/qqiu/project/others/test/integration/cm_VR.snp_gene.out", sep = "\t", header = T)

mc_NV_dis <- mc_NV %>%
  dplyr::select(SNP, gene_id, gene_name, distance_category) %>%
  mutate(distance_category = factor(distance_category, 
                                    levels=c("<10 kb", "10-50 kb", "50-100 kb", "100-200 kb", "200-500 kb", "> 500kb")),
         sample = "Arteriole", method = "Micro-C") %>%
  unique()
mc_VR_dis <- mc_VR %>%
  dplyr::select(SNP, gene_id, gene_name, distance_category) %>%
  mutate(distance_category = factor(distance_category, 
                                    levels=c("<10 kb", "10-50 kb", "50-100 kb", "100-200 kb", "200-500 kb", "> 500kb")),
         sample = "EDA", method = "Micro-C") %>%
  unique()
cm_NV_dis <- cm_NV %>%
  dplyr::select(SNP, gene_id, gene_name, distance_category) %>%
  mutate(distance_category = factor(distance_category, 
                                    levels=c("<10 kb", "10-50 kb", "50-100 kb", "100-200 kb", "200-500 kb", "> 500kb")),
         sample = "Arteriole", method = "Capture Micro-C") %>%
  unique()
cm_VR_dis <- cm_VR %>%
  dplyr::select(SNP, gene_id, gene_name, distance_category) %>%
  mutate(distance_category = factor(distance_category, 
                                    levels=c("<10 kb", "10-50 kb", "50-100 kb", "100-200 kb", "200-500 kb", "> 500kb")),
         sample = "EDA", method = "Capture Micro-C") %>%
  unique()

snp_gene_dis = rbind(mc_NV_dis, mc_VR_dis, cm_NV_dis, cm_VR_dis)
snp_gene_dis$method = factor(snp_gene_dis$method, levels=c("Micro-C", "Capture Micro-C"))
p = ggplot(snp_gene_dis, aes(x = distance_category, fill = sample)) +
  geom_bar(width=.8, position = "dodge") +
  scale_fill_manual(values = c("Arteriole"="skyblue2", "EDA"="salmon2")) +
  theme_classic() +
  theme(strip.text = element_text(size = 10),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "SNP distance to gene TSS",
       y = "Number of SNP-gene pairs",
       fill = "Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~method, scales = "free")
p
ggsave("/xdisk/mliang1/qqiu/project/others/test/figure/fig4d.snp_gene_dis.barplot.png", width=482/96, height=242/96, dpi=300)




### functional enrichment of genes interacted with SNPs
expr <- read.csv("/xdisk/mliang1/qqiu/project/others/test/data/vessel.gene.FPKM.rm_013.csv", header = T, row.names = 1)

input_df <- data.frame(input = c("mc_NV", "mc_VR", "cm_NV", "cm_VR"),
                      title = c("Micro-C (arteriole)", "Micro-C (EDA)", "Capture Micro-C (arteriole)", "Capture Micro-C (EDA)"))
for( i in 1:nrow(input_df)){
  
  dataset <- input_df$input[i]
  title <- input_df$title[i]
  input <- paste0("/xdisk/mliang1/qqiu/project/others/test/integration/", dataset, ".snp_gene.out")
  outfile <- paste0("/xdisk/mliang1/qqiu/project/others/test/figure/fig4e.", dataset, ".snp_gene.GO_enrich.png")
  
  if(grepl("NV", dataset)){
    expr_list <- colnames(expr)[grep("namR$", colnames(expr))]
  }else{
    expr_list <- colnames(expr)[grep("vrmR$", colnames(expr))]
  }
  expr_gene_list <- rownames(expr[rowMeans(expr[expr_list])>0.1,])
  
  snp_gene <- read.table(input, header = T, sep = "\t")
  gene_list <- unique(snp_gene$gene_id)
  gene_entrez <- bitr(gene_list, fromType = "ENSEMBL", 
                      toType = "ENTREZID", 
                      OrgDb = org.Hs.eg.db)
  
  go_enrichment <- enrichGO(gene = gene_entrez$ENTREZID,
                            OrgDb = org.Hs.eg.db,
                            # universe = gene_entrez$ENTREZID,
                            keyType = "ENTREZID",
                            ont = "BP", 
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05)
  
  if(nrow(go_enrichment)>0){
    p = barplot(go_enrichment, showCategory = 10, x = "GeneRatio", title = title)
    print(p)
    ggsave(outfile, width=546/96, height=493/96, dpi=300)
  }
  
}





### SNP-gene interact pattern 
perform_fisher_test <- function(Freq_obs, TYPE, Bin, data) {
  
  cell_value <- Freq_obs
  same_type_uncovered <- sum(data[data$TYPE == TYPE & data$Bin == Bin, ]$Uncovered_SNPs)
  diff_type_covered <- sum(data[data$TYPE != TYPE & data$Bin == Bin, ]$Covered_SNPs)
  diff_type_uncovered <- sum(data[data$TYPE != TYPE & data$Bin == Bin, ]$Uncovered_SNPs)
  
  contingency_table <- matrix(c(cell_value, same_type_uncovered, 
                                diff_type_covered, diff_type_uncovered),
                              nrow = 2, byrow = TRUE)
  
  if (any(contingency_table < 0) || sum(contingency_table) == 0 || any(is.na(contingency_table))) {
    return(list(p_value = NA, odds_ratio = NA))
  }
  
  fisher_result <- tryCatch({
    fisher.test(contingency_table)
  }, error = function(e) {
    return(list(p.value = NA, estimate = NA))
  })
  
  return(list(p_value = fisher_result$p.value, odds_ratio = fisher_result$estimate))
}



trait_list <- c("systolic_bp", "diastolic_bp", "pulse_pressure", "essential_hypertension", "stroke", 
               "diabetic_nephropathy", "diabetic_retinopathy", "peripheral_arterial_disease")
gwas <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/gwas_catalog.processed.txt", header = T, sep = "\t")
gwas <- gwas[gwas$TRAIT %in% trait_list, ]
gwas[is.na(gwas$REG_ANNO), ]$REG_ANNO <- "unknown"
gwas_reg <- unique(gwas[, c("SNP", "REG_ANNO")])

input_df = data.frame(input = c("cm_NV", "cm_VR"),
                      title = c("Arteriole", "EDA"))

for( i in 1:nrow(input_df)){
  
  dataset <- input_df$input[i]
  input <- paste0("/xdisk/mliang1/qqiu/project/others/test/integration/", dataset, ".snp_gene.out")
  snp_gene <- read.table(input, header = T, sep = "\t")
  outfile <- paste0(input_df[i,]$input, ".snp_by_reg.bar_plot.png")
  title <- input_df[i,]$title
  
  total_snps <- data.frame(TYPE = names(table(gwas_reg$REG_ANNO)),
                           Total_SNPs = as.numeric(table(gwas_reg$REG_ANNO)))
  
  BIN1_SNP <- snp_gene[snp_gene$BIN_ID_SNP=="BIN1",]$SNP
  covered_snps_bin1 <- data.frame(TYPE = names(table(gwas_reg[gwas_reg$SNP %in% BIN1_SNP, ]$REG_ANNO)),
                             Covered_SNPs = as.numeric(table(gwas_reg[gwas_reg$SNP %in% BIN1_SNP, ]$REG_ANNO)))
  BIN2_SNP <- snp_gene[snp_gene$BIN_ID_SNP=="BIN2",]$SNP
  covered_snps_bin2 <- data.frame(TYPE = names(table(gwas_reg[gwas_reg$SNP %in% BIN2_SNP, ]$REG_ANNO)),
                             Covered_SNPs = as.numeric(table(gwas_reg[gwas_reg$SNP %in% BIN2_SNP, ]$REG_ANNO)))

  covered_rate_bin1 <- length(unique(intersect(gwas_reg$SNP, BIN1_SNP)))/length(unique(gwas_reg$SNP))
  covered_rate_bin2 <- length(unique(intersect(gwas_reg$SNP, BIN2_SNP)))/length(unique(gwas_reg$SNP))
  
  covered_rate <- data.frame(Bin = c("Promoter bin", "Distal bin"),
                            rate = c(covered_rate_bin1, covered_rate_bin2))
  
  merged_data_bin1 <- merge(total_snps, covered_snps_bin1, by = "TYPE", all.x = TRUE) %>%
    mutate(Covered_SNPs = ifelse(is.na(Covered_SNPs), 0, Covered_SNPs),
           Uncovered_SNPs = Total_SNPs - Covered_SNPs,
           Bin = "Promoter bin") %>%
    mutate(Proportion = Covered_SNPs / Total_SNPs)

  merged_data_bin2 <- merge(total_snps, covered_snps_bin2, by = "TYPE", all.x = TRUE) %>%
    mutate(Covered_SNPs = ifelse(is.na(Covered_SNPs), 0, Covered_SNPs),
           Uncovered_SNPs = Total_SNPs - Covered_SNPs,
           Bin = "Distal bin") %>%
    mutate(Proportion = Covered_SNPs / Total_SNPs)
  
  merged_data <- rbind(merged_data_bin1, merged_data_bin2)
  merged_data <- merged_data %>%
    rowwise() %>%
    mutate(result = list(perform_fisher_test(Covered_SNPs, TYPE, Bin, merged_data))) %>%
    unnest_wider(result) %>%
    ungroup() %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    mutate(
      significance = case_when(
        p_adj < 0.001 & odds_ratio>1 ~ "***",
        p_adj < 0.01 & odds_ratio>1 ~ "**",
        p_adj < 0.05 & odds_ratio>1 ~ "*",
        TRUE ~ ""
      )
    )
  
  p1 <- ggplot(merged_data, aes(x = reorder(TYPE, Total_SNPs), y = Proportion)) +
    geom_bar(stat = "identity", fill="steelblue") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    geom_hline(data = covered_rate, aes(yintercept = rate), linetype = "dashed", color = "red", size = 1) +
    geom_text(aes(label = significance), # y = Proportion+0.1*max(Proportion)), 
              color = "black", size = 5) +
    labs(title = title,
         x = "SNP regulatory annotation",
         y = "Proportion of covered SNPs",
         fill = "SNP Status") +
    theme_classic() +
    theme(axis.text = element_text(colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, size=10),
          legend.position = "none") +
    facet_grid(fct_rev(Bin)~.)
  
  print(p1)
  
  ggsave(paste0("/xdisk/mliang1/qqiu/project/others/test/figure/fig4g.", outfile), width=282/96, height=364/96, dpi=300)
  
}





### single-cell rna-seq

cm_NV <- read.table("/xdisk/mliang1/qqiu/project/others/test/integration/cm_NV.snp_gene.out", header = T, sep = "\t")
cm_VR <- read.table("/xdisk/mliang1/qqiu/project/others/test/integration/cm_VR.snp_gene.out", header = T, sep = "\t")

cm_NV[is.na(cm_NV$REG_ANNO), ]$REG_ANNO <- "unknown"
cm_NV <- unique(cm_NV[, c("SNP", "gene_id", "gene_name", "distance_category", "is_TF", "REG_ANNO", "overlap_type", "mc_id", "score")])
cm_VR[is.na(cm_VR$REG_ANNO), ]$REG_ANNO <- "unknown"
cm_VR <- unique(cm_VR[, c("SNP", "gene_id", "gene_name", "distance_category", "is_TF", "REG_ANNO", "overlap_type", "mc_id", "score")])

merged_data <- merge(cm_NV, cm_VR, by = c("SNP", "gene_id", "gene_name", "distance_category", "is_TF", "REG_ANNO"), suffixes = c("_NV", "_VR"), all=T)
merged_data <- merged_data %>%
  mutate(score_NV = ifelse(is.na(score_NV), 0, score_NV),
         score_VR = ifelse(is.na(score_VR), 0, score_VR),
         score_diff = score_NV - score_VR)


filtered_data <- merged_data
overall_gene_list <- unique(filtered_data$gene_name)
filtered_data <- filtered_data %>%
  filter(score_NV*score_VR==0) %>%
  mutate(specificity = ifelse(score_NV==0, "EDA", "Arteriole"))

arteriole_gene_list <- unique(filtered_data[filtered_data$specificity=="Arteriole",]$gene_name)
eda_gene_list <- unique(filtered_data[filtered_data$specificity=="EDA",]$gene_name)

arteriole_specific_gene_list <- setdiff(arteriole_gene_list, eda_gene_list)
eda_specific_gene_list <- setdiff(eda_gene_list, arteriole_gene_list)
overlap_gene_list <- setdiff(overall_gene_list, c(arteriole_specific_gene_list, eda_specific_gene_list))

full_specific_gene_list <- data.frame(gene = c(arteriole_specific_gene_list, eda_specific_gene_list, overlap_gene_list),
                                specificity = c(rep("Arteriole specific", length(arteriole_specific_gene_list)), 
                                                rep("EDA specific", length(eda_specific_gene_list)),
                                                rep("Overlap", length(overlap_gene_list))))

filtered_data <- merged_data %>%
  filter(REG_ANNO=="enhancer", (overlap_type_NV %in% c("promoter") | 
                                  overlap_type_VR %in% c("promoter") )) 
overall_gene_list <- unique(filtered_data$gene_name)
filtered_data <- filtered_data%>%
  filter(score_NV*score_VR==0) %>%
  mutate(specificity = ifelse(score_NV==0, "EDA", "Arteriole"))

arteriole_gene_list <- unique(filtered_data[filtered_data$specificity=="Arteriole",]$gene_name)
eda_gene_list <- unique(filtered_data[filtered_data$specificity=="EDA",]$gene_name)

arteriole_specific_gene_list <- setdiff(arteriole_gene_list, eda_gene_list)
eda_specific_gene_list <- setdiff(eda_gene_list, arteriole_gene_list)
overlap_gene_list <- setdiff(overall_gene_list, c(arteriole_specific_gene_list, eda_specific_gene_list))

espg_specific_gene_list <- data.frame(gene = c(arteriole_specific_gene_list, eda_specific_gene_list, overlap_gene_list),
                                     specificity = c(rep("Arteriole specific", length(arteriole_specific_gene_list)), 
                                                     rep("EDA specific", length(eda_specific_gene_list)),
                                                     rep("Overlap", length(overlap_gene_list))))


filtered_data <- merged_data %>%
  filter(REG_ANNO=="CTCF_binding_site", (overlap_type_NV %in% c("promoter") | 
                                  overlap_type_VR %in% c("promoter") )) 
overall_gene_list <- unique(filtered_data$gene_name)
filtered_data <- filtered_data %>%
    filter(score_NV*score_VR==0) %>%
  mutate(specificity = ifelse(score_NV==0, "EDA", "Arteriole"))

arteriole_gene_list <- unique(filtered_data[filtered_data$specificity=="Arteriole",]$gene_name)
eda_gene_list <- unique(filtered_data[filtered_data$specificity=="EDA",]$gene_name)

arteriole_specific_gene_list <- setdiff(arteriole_gene_list, eda_gene_list)
eda_specific_gene_list <- setdiff(eda_gene_list, arteriole_gene_list)
overlap_gene_list <- setdiff(overall_gene_list, c(arteriole_specific_gene_list, eda_specific_gene_list))

CTCFpg_specific_gene_list <- data.frame(gene = c(arteriole_specific_gene_list, eda_specific_gene_list, overlap_gene_list),
                                       specificity = c(rep("Arteriole specific", length(arteriole_specific_gene_list)), 
                                                       rep("EDA specific", length(eda_specific_gene_list)),
                                                       rep("Overlap", length(overlap_gene_list))))


filtered_data <- merged_data %>%
  filter((overlap_type_NV %in% c("promoter") | 
            overlap_type_VR %in% c("promoter") )) 
overall_gene_list <- unique(filtered_data$gene_name)
filtered_data <- filtered_data %>%
  filter(score_NV*score_VR==0) %>%
  mutate(specificity = ifelse(score_NV==0, "EDA", "Arteriole"))

arteriole_gene_list <- unique(filtered_data[filtered_data$specificity=="Arteriole",]$gene_name)
eda_gene_list <- unique(filtered_data[filtered_data$specificity=="EDA",]$gene_name)

arteriole_specific_gene_list <- setdiff(arteriole_gene_list, eda_gene_list)
eda_specific_gene_list <- setdiff(eda_gene_list, arteriole_gene_list)
overlap_gene_list <- setdiff(overall_gene_list, c(arteriole_specific_gene_list, eda_specific_gene_list))

pg_specific_gene_list <- data.frame(gene = c(arteriole_specific_gene_list, eda_specific_gene_list, overlap_gene_list),
                                   specificity = c(rep("Arteriole specific", length(arteriole_specific_gene_list)), 
                                                   rep("EDA specific", length(eda_specific_gene_list)),
                                                   rep("Overlap", length(overlap_gene_list))))



seurat_object = readRDS("/xdisk/mliang1/qqiu/project/others/test/cluster/huves1SN.cluster.rds")
ta_obj = readRDS("/xdisk/mliang1/qqiu/project/others/test/data/NG-2022/ascending_descending_human_aorta.rds")

# arteriole_avg = AverageExpression(seurat_object, return.seurat = F, add.ident = 'cell_type')$RNA
# ta_avg = AverageExpression(seurat_object, return.seurat = F, add.ident = 'cell_type')$RNA

# arteriole_markers = FindAllMarkers(seurat_object, group.by="cell_type", assay = "RNA", only.pos = TRUE, min.pct = 0.1)
# arteriole_markers$pct.diff = arteriole_markers$pct.1 - arteriole_markers$pct.2
# arteriole_markers$sample = "Arteriole"
# 
# ta_markers = FindAllMarkers(ta_obj, group.by="cell_type", assay = "RNA", only.pos = TRUE, min.pct = 0.1)
# ta_markers$pct.diff = ta_markers$pct.1 - ta_markers$pct.2
# ta_markers$sample = "Thoracic aorta"
# save(arteriole_markers, ta_markers, file = "/xdisk/mliang1/qqiu/project/others/test/data/snRNAseq.deg.RData")


load("/xdisk/mliang1/qqiu/project/others/test/data/snRNAseq.deg.RData")

markers_merged = arteriole_markers %>%
  filter(p_val_adj<0.05 & avg_log2FC>0.25)

total_genes = rownames(seurat_object)


fisher_test_function <- function(specific_gene_list, markers_merged, total_genes){
  
  results <- data.frame()
  for (capture_spec in unique(specific_gene_list$specificity)) {
    
    capture_genes_subset <- specific_gene_list %>%
      filter(specificity == capture_spec) %>%
      pull(gene)
    
    cat("\nCapture specificity:", capture_spec, "- Number of genes:", length(capture_genes_subset), "\n")
    
    for (si in unique(markers_merged$sample)) {
      for (ci in unique(markers_merged$cluster)) {
        
        markers_subset <- markers_merged %>%
          filter(sample == si, cluster == ci) %>%
          pull(gene)
        
        overlap_genes <- intersect(capture_genes_subset, markers_subset)
        
        a <- length(overlap_genes)  
        b <- length(capture_genes_subset) - a  
        c <- length(markers_subset) - a 
        d <- length(total_genes) - (a + b + c) 
        
        fisher_matrix <- matrix(c(a, b, c, d), nrow = 2)
        fisher_result <- fisher.test(fisher_matrix)
        
        odds_ratio <- fisher_result$estimate
        p_value <- fisher_result$p.value
        
        results <- rbind(results, data.frame(
          capture_microC_specificity = capture_spec,
          sample = si,
          cell_type = ci,
          odds_ratio = odds_ratio,
          p_value = p_value
        ))
      }
    }
  }
  
  results$p_adj <- p.adjust(results$p_value, "BH")
  return(results)
}

full_fisher_test_res <- fisher_test_function(full_specific_gene_list, markers_merged, total_genes)
espg_fisher_test_res <- fisher_test_function(espg_specific_gene_list, markers_merged, total_genes)
pg_fisher_test_res <- fisher_test_function(pg_specific_gene_list, markers_merged, total_genes)

full_fisher_test_res$gene_type <- "All interacting genes"
pg_fisher_test_res$gene_type <- "Promoter-interacting genes"
espg_fisher_test_res$gene_type <- "Promoter-enhancer SNP\ninteracting genes"

fisher_test_res <- rbind(full_fisher_test_res,
                         pg_fisher_test_res,
                         espg_fisher_test_res) %>%
  mutate(gene_type = factor(gene_type, levels=c("Promoter-enhancer SNP\ninteracting genes",
                                                "Promoter-interacting genes",
                                                "All interacting genes"))) %>%
  mutate(capture_microC_specificity = factor(capture_microC_specificity, 
                                             levels=c("Overlap",
                                                      "EDA specific",
                                                      "Arteriole specific"))) %>%
  filter(odds_ratio != 0) %>%
  mutate(log_p_value = -log10(p_value),
         significant = ifelse(p_value < 0.05, "Significant", "Not Significant"))

p <- ggplot(fisher_test_res, aes(x = cell_type, 
                                 y = capture_microC_specificity, 
                                 size = log_p_value, 
                                 color = odds_ratio)) +
  geom_point() + 
  geom_point(data = fisher_test_res %>% filter(significant == "Significant"), 
             aes(x = cell_type, 
                 y = capture_microC_specificity), 
             shape = 21, fill = NA, color = "black", stroke = 1) + 
  theme_classic() +
  scale_color_gradient2(low = "blue", high = "red", mid = "white", midpoint = 1) +  # Adjust color gradient for odds ratio
  labs(
    x = "snRNA-seq cell type",
    y = "Capture Micro-C\nspecificity",
    size = "-log10(p-value)",
    color = "Odds Ratio",
    subtitle = "Enrichment of interacting genes in cell type specific genes"
  ) +
  theme(
    strip.text = element_text(color="black", size = 10),
    axis.text = element_text(color="black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    legend.box = "horizontal"
  ) +
  facet_wrap(~gene_type, scales = "free_x")
print(p)

ggsave("/xdisk/mliang1/qqiu/project/others/test/figure/fig4h.interacting_gene_cell_type_specificity.png", width=896/96, height=219/96, dpi=300)




fib_arteriole <- fib_arteriole <- markers_merged %>%
  filter(sample == "Arteriole", cluster == "Fibroblast") %>%
  left_join(merged_data %>% 
              filter(REG_ANNO == "enhancer", overlap_type_NV == "promoter"), 
            by = c("gene" = "gene_name")) %>%
  mutate(
    specificity = case_when(
      score_NV != 0 & score_VR == 0 ~ "Arteriole specific",
      score_NV == 0 & score_VR != 0 ~ "EDA specific", 
      score_NV != 0 & score_VR != 0 ~ "Overlap" 
    )
  ) %>%
  arrange(desc(pct.diff)) %>%
  mutate(rank = row_number())

input_dat <- fib_arteriole

ranked_genes <- input_dat %>%
  arrange(desc(avg_log2FC)) %>%
  mutate(rank = row_number())

highlighted_genes <- ranked_genes %>%
  filter(!is.na(SNP) & !is.na(gene_id)) %>%
  group_by(gene) %>% 
  slice(1) %>%
  ungroup() %>%
  arrange(desc(avg_log2FC)) %>%
  mutate(snp_gene_label = paste(SNP, gene, sep = "-"))

p <- ggplot(ranked_genes, aes(x = rank, y = avg_log2FC)) +
  geom_point(color = "grey") +
  geom_point(data = highlighted_genes, aes(x = rank, y = avg_log2FC, color = specificity), size = 3) +
  geom_text_repel(data = highlighted_genes[1:10,], aes(x = rank, y = avg_log2FC, label = snp_gene_label),
                  size = 4, color = "black", nudge_y = c(-1,1), nudge_y = c(-1,1), box.padding = 0.5, 
                  point.padding = 0.5, segment.color = "grey50", max.overlaps = Inf,
                  force = 20) +
  scale_color_manual(values = c("skyblue2", "black")) +
  labs(title = "Fibroblast specific genes",
       x = "Rank", color = "Capture Micro-C\nspecificity",
       y = "Log2 fold-change (compared to other cells)") +
  theme_classic() +
  theme(axis.text = element_text(color = "black"))

print(p)
ggsave("/xdisk/mliang1/qqiu/project/others/test/figure/fig4i.fibroblast_gene.snp.png", width=622/96, height=333/96, dpi=300)


gene <- c("MEG3", "SEMA4A")
p1 <- seurat_object %>% FeaturePlot(., features = gene, order = T) & labs(x="", y="") & 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), plot.margin = unit(c(0.1, 0, 0, 0), "cm"))
print(p1)
ggsave("/xdisk/mliang1/qqiu/project/others/test/figure/fig4j.fibroblast.selected_gene.umap.png", width=384/96, height=165/96, dpi=300)






################################################################################
### SNP-gene nomination
### preprocess snp-gene pairs by cap micro-c
library(DiagrammeR)


cm_NV <- read.table("/xdisk/mliang1/qqiu/project/others/test/integration/cm_NV.snp_gene.out", header = T, sep = "\t")
cm_VR <- read.table("/xdisk/mliang1/qqiu/project/others/test/integration/cm_VR.snp_gene.out", header = T, sep = "\t")

cm_NV[is.na(cm_NV$REG_ANNO), ]$REG_ANNO <- "unknown"
cm_NV <- unique(cm_NV[, c("SNP", "SNP_CHR", "SNP_POS", "gene_id", "gene_name", "distance_category", "is_TF", "REG_ANNO", "overlap_type", "mc_id", "score")])
cm_VR[is.na(cm_VR$REG_ANNO), ]$REG_ANNO <- "unknown"
cm_VR <- unique(cm_VR[, c("SNP", "SNP_CHR", "SNP_POS", "gene_id", "gene_name", "distance_category", "is_TF", "REG_ANNO", "overlap_type", "mc_id", "score")])

snp_gene_merged <- rbind(cm_NV, cm_VR)

### conservation
hg38_mm10 <- rtracklayer::import.chain("/xdisk/mliang1/qqiu/reference/liftover/hg38ToMm10.over.chain")
hg38_rn7 <- rtracklayer::import.chain("/xdisk/mliang1/qqiu/reference/liftover/hg38ToRn7.over.chain")

SNPS <- snp_gene_merged %>% 
  dplyr::select(SNP, SNP_CHR, SNP_POS) %>% unique()

gr <- GRanges(seqnames = SNPS$SNP_CHR,
             ranges = IRanges(start = SNPS$SNP_POS,
                              end = SNPS$SNP_POS),
             POS = SNPS$SNP_POS)
snp_mm10 <- rtracklayer::liftOver(x = gr, chain = hg38_mm10) %>% unlist()
snp_rn7 <- rtracklayer::liftOver(x = gr, chain = hg38_rn7) %>% unlist()

lifted_mm10 <- data.frame(original_POS = mcols(snp_mm10)$POS,
                          POS_mm10 = paste0(seqnames(snp_mm10), "_", 
                                            start(snp_mm10), "_", 
                                            end(snp_mm10)))

lifted_rn7 <- data.frame(original_POS = mcols(snp_rn7)$POS,
                         POS_mm10 = paste0(seqnames(snp_rn7), "_", 
                                           start(snp_rn7), "_", 
                                           end(snp_rn7)))


### LD info
snps_LD_1 <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/SNPS_in_LD.txt", header = T, sep = "\t")
snps_LD_2 <- read.csv("/xdisk/mliang1/qqiu/project/others/test/data/LD_info_Joan_092724.csv", header = T)

ld_1 <- c(names(table(snps_LD_1$rsID1)[table(snps_LD_1$rsID1)==1]), snps_LD_2[snps_LD_2$total_SNPs_inLD==1, ]$sentinal_SNP) %>% unique()


### BP trait list
gwas <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/gwas_catalog.processed.txt", header = T, sep = "\t")
bp_snp_list <- unique(gwas[gwas$TRAIT %in% c("diastolic_bp", "pulse_pressure", "systolic_bp"),]$SNP)

# bp_snp_list = read.table("/xdisk/mliang1/qqiu/data/HT-GWAS/BP_SNP.NG-2024.out", header = T, sep = "\t")
# bp_snp_list = bp_snp_list$rsID

### protein coding list
ensembl <- read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.GRCH38.out", sep = '\t', header=T)
pc_gene_list <- unique(ensembl[ensembl$Gene.type %in% c("protein_coding"),]$Gene.stable.ID)


### low meth list
meth <- read.csv("/xdisk/mliang1/qqiu/project/others/test/data/vessel.RRBS.dat.cln.promoter.methyl.out", header = T, row.names = 1)
meth_gene_list <- rownames(meth[rowMeans(meth, na.rm = T)<0.2, ])


### expr gene list
expr = read.csv("/xdisk/mliang1/qqiu/project/others/test/data/vessel.gene.FPKM.rm_013.csv", header = T, row.names = 1)
expr_gene_list <- rownames(expr[rowMeans(log2(expr+1))>1,])


### filtration
snp_gene_df_filter <- snp_gene_merged %>% 
  dplyr::select(SNP, SNP_CHR, SNP_POS, gene_id, gene_name, REG_ANNO, overlap_type, distance_category) %>% 
  group_by(across(-c(REG_ANNO, overlap_type))) %>%
  summarise(overlap_type = paste(sort(unique(overlap_type)), collapse = ", "),
            REG_ANNO = paste(sort(unique(REG_ANNO)), collapse = ", ")) %>% 
  unique() %>%
  filter(SNP %in% ld_1, 
         SNP_POS %in% lifted_mm10$original_POS &
           SNP_POS %in% lifted_rn7$original_POS,
         # trait %in% c("diastolic_bp", "systolic_bp", "pulse_pressure"),
         distance_category %in% c("100-200 kb", "200-500 kb", "> 500kb")) %>%
  group_by(gene_id) %>% dplyr::mutate(n_SNP_per_gene = n_distinct(SNP)) %>% ungroup() %>%
  filter(n_SNP_per_gene==1,
         grepl("promoter", overlap_type),
         SNP %in% bp_snp_list,
         gene_id %in% pc_gene_list)




filter_counts <- list()
snp_gene_df_filter <- snp_gene_merged %>% 
  dplyr::select(SNP, SNP_CHR, SNP_POS, gene_id, gene_name, REG_ANNO, overlap_type, distance_category, mc_id, score) %>% 
  group_by(across(-c(REG_ANNO, overlap_type, mc_id, score))) %>%
  summarise(overlap_type = paste(sort(unique(overlap_type)), collapse = ", "),
            REG_ANNO = paste(sort(unique(REG_ANNO)), collapse = ", "),
            mc_id = paste(sort(unique(mc_id)), collapse = ", "),
            score = paste(sort(unique(score)), collapse = ", ")) %>% 
  unique()
snp_gene_df_filter <- snp_gene_df_filter %>%
  filter(grepl("promoter", overlap_type))
filter_counts[["Initial"]] <- nrow(snp_gene_df_filter)

snp_gene_df_filter <- snp_gene_df_filter %>%
  filter(SNP %in% bp_snp_list)
filter_counts[["BP associated"]] <- nrow(snp_gene_df_filter)

snp_gene_df_filter <- snp_gene_df_filter %>%
  filter(gene_id %in% pc_gene_list)
filter_counts[["Interact with promoter of protein-coding genes"]] <- nrow(snp_gene_df_filter)

snp_gene_df_filter <- snp_gene_df_filter %>%
  filter(gene_id %in% meth_gene_list,
         gene_id %in% expr_gene_list)
filter_counts[["Low promoter methylation and\nminimum gene expression"]] <- nrow(snp_gene_df_filter)

snp_gene_df_filter <- snp_gene_df_filter %>%
  filter(SNP %in% ld_1)
filter_counts[["SNPs isolated within LD blocks"]] <- nrow(snp_gene_df_filter)

snp_gene_df_filter <- snp_gene_df_filter %>%
  filter(SNP_POS %in% lifted_mm10$original_POS &
           SNP_POS %in% lifted_rn7$original_POS)
filter_counts[["Cross-species conservation"]] <- nrow(snp_gene_df_filter)

snp_gene_df_filter <- snp_gene_df_filter %>%
  group_by(gene_id) %>%
  dplyr::mutate(n_SNP_per_gene = n_distinct(SNP)) %>%
  ungroup() %>%
  filter(n_SNP_per_gene == 1)
filter_counts[["Genes paired with only one SNP"]] <- nrow(snp_gene_df_filter)

snp_gene_df_filter <- snp_gene_df_filter %>%
  filter(distance_category %in% c("100-200 kb", "200-500 kb", "> 500kb"))
filter_counts[["Distal regulatory element (>100 kb)"]] <- nrow(snp_gene_df_filter)

filter_data <- data.frame(Step = names(filter_counts),
                         SNP_Gene_Count = as.numeric(filter_counts))

generate_grViz <- function(filter_data, ranksep = 0.25, nodesep = 0.5) {
  nodes <- paste0(filter_data$Step, ":\n", filter_data$SNP_Gene_Count, " SNP-gene pairs")
  edges <- paste0(seq(1, nrow(filter_data) - 1), " -> ", seq(2, nrow(filter_data)))
  
  diagram_code <- paste0(
    "digraph G {
       graph [layout = dot, rankdir = TB, ranksep = ", ranksep, ", nodesep = ", nodesep, "]
       node [shape = box, style = filled, fillcolor = lightblue, fontname = Helvetica]
       ",
    paste0(paste0(seq(1, nrow(filter_data)), " [label = \"", nodes, "\"]"), collapse = "\n"),
    "\n",
    paste0(edges, collapse = "\n"),
    "\n}"
  )
  grViz(diagram_code)
}

generate_grViz(filter_data)

