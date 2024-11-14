library(dplyr)
library(stringr)
library(ggVennDiagram)
library(grid)
library(GenomicRanges)


################################################################################
### functions
standardize_column_names <- function(df) {
  df %>%
    rename_with(~"BIN1_CHR",  matches("bait_chr|bin1_chr|chr1|BIN1_CHR")) %>%
    rename_with(~"BIN1_START", matches("bait_start|bin1_start|start1|BIN1_START")) %>%
    rename_with(~"BIN1_END", matches("bait_end|bin1_end|end1|BIN1_END")) %>%
    
    rename_with(~"BIN2_CHR",  matches("otherEnd_chr|bin2_chr|chr2|BIN2_CHR|BIN2_CHROMOSOME")) %>%
    rename_with(~"BIN2_START", matches("otherEnd_start|bin2_start|start2|BIN2_START")) %>%
    rename_with(~"BIN2_END", matches("otherEnd_end|bin2_end|end2|BIN2_END"))
}



################################################################################
### preprocessing original probes
ref <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/Homo_sapiens.GRCh38.112.transcript.gtf", header = F, sep = '\t')
ref$transcript_id <- str_extract(ref$V9, "transcript_id [^;]+") %>%
  gsub("transcript_id ", "", .)
ref$gene_id <- str_extract(ref$V9, "gene_id [^;]+") %>%
  gsub("gene_id ", "", .)

ref_ranges <- GRanges(
  seqnames = Rle(paste0("chr", ref$V1)),
  ranges = IRanges(start = ref$V4, end = ref$V5),
  strand = Rle(ref$V7),
  transcript_id = ref$transcript_id,
  gene_id = ref$gene_id
)
promoter_ranges <- promoters(ref_ranges, upstream = 1000, downstream = 500)

probes <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/h_probes_v1.0.bed", header = F)

probe_ranges <- GRanges(
  seqnames = probes$V1,
  ranges = IRanges(start = probes$V2, end = probes$V3)
)

probe_hits <- findOverlaps(probe_ranges, promoter_ranges)
probe_covered <- promoter_ranges[subjectHits(probe_hits), ]

write.table(probe_covered, "/xdisk/mliang1/qqiu/project/others/test/data/h_probes_v1.0.trans_gene.out", sep = ",", row.names = F, quote = F)



################################################################################
### load expr file and process
library(DESeq2)
count_data <- read.csv("/xdisk/mliang1/qqiu/project/others/test/data/vessel.gene_count_matrix.cln.csv", header = TRUE, row.names = 1, sep = ",")
count_data <- count_data[, !(grepl("013", colnames(count_data)))]
rownames(count_data) <- gsub("\\..*", "", sapply(rownames(count_data), function(x) strsplit(x, "\\|")[[1]][1]))
condition <- factor(rep(c("EC", "NA", "VR"), 3))

dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = data.frame(condition),
                              design = ~ 1)

dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

write.csv(normalized_counts, file = "/xdisk/mliang1/qqiu/project/others/test/data/vessel.gene.normalized_counts.rm_013.csv")


FPKM_data <- read.csv("/xdisk/mliang1/qqiu/project/others/test/data/FPKM.txt.cln.tsv", header = TRUE, sep = "\t")
FPKM_data <- FPKM_data[, !(grepl("013", colnames(FPKM_data)))]
FPKM_data$Gene_ID <- gsub("\\..*", "", FPKM_data$Gene_ID)
FPKM_data_agg <- FPKM_data %>%
  group_by(Gene_ID) %>%
  summarise(across(everything(), mean)) %>%
  as.data.frame()
rownames(FPKM_data_agg) <- FPKM_data_agg$Gene_ID
FPKM_data_agg <- FPKM_data_agg[, grepl("mR", colnames(FPKM_data_agg))]
write.csv(FPKM_data_agg, file = "/xdisk/mliang1/qqiu/project/others/test/data/vessel.gene.FPKM.rm_013.csv")


################################################################################
### preprocessing interaction files
filter_promoter_regions <- function(input_file, promoter_ranges) {
  
  output_file = gsub("csv", "filtered.csv", input_file)
  
  input_df <- fread(input_file)
  
  cm_ranges <- GRanges(
    seqnames = input_df$bait_chr,
    ranges = IRanges(start = input_df$bait_start, end = input_df$bait_end)
  )
  
  bait_hits <- findOverlaps(cm_ranges, promoter_ranges)
  unique_hits <- unique(queryHits(bait_hits))
  filtered_df <- input_df[unique_hits, ]
  
  write.table(filtered_df, file = output_file, sep = ",", row.names = F)
  
  cat("Original records:", nrow(input_df), "\n",
      "Filtered records:", nrow(filtered_df), "\n",
      "Convert rate:", nrow(filtered_df)/nrow(input_df), "\n")
}


ref <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/Homo_sapiens.GRCh38.112.transcript.gtf", header = F, sep = '\t')
ref_ranges <- GRanges(
  seqnames = Rle(paste0("chr", ref$V1)),
  ranges = IRanges(start = ref$V4, end = ref$V5),
  strand = Rle(ref$V7)
)
promoter_ranges <- promoters(ref_ranges, upstream = 1000, downstream = 500)

filter_promoter_regions("/xdisk/mliang1/qqiu/project/others/test/data/NVCMC_newpara_10kb_intra.csv", promoter_ranges)
filter_promoter_regions("/xdisk/mliang1/qqiu/project/others/test/data/NVCMC_newpara_20kb_intra.csv", promoter_ranges)
filter_promoter_regions("/xdisk/mliang1/qqiu/project/others/test/data/VRCMC_newpara_10kb_intra.csv", promoter_ranges)
filter_promoter_regions("/xdisk/mliang1/qqiu/project/others/test/data/VRCMC_newpara_20kb_intra.csv", promoter_ranges)



################################################################################
### venn plot for loops and interactions
venn_plot <- function(df1, df2) {
  
  sample1 <- ifelse(grepl("NV", df1), "Arterioles", "EDA")
  sample2 <- ifelse(grepl("NV", df2), "Arterioles", "EDA")
  
  inter <- ifelse(grepl("^mc", df1), "Loops", "Interactions")
  resolution <- strsplit(df1, "_")[[1]][3] %>% gsub("kb", " kb", .)
  plot_title <- paste0(inter, " (", resolution, ")\n")
  
  df1 <- standardize_column_names(get(df1))
  df2 <- standardize_column_names(get(df2))
  
  df1_count <- nrow(df1)
  df2_count <- nrow(df2)
  
  overlap_count <- 0
  overlap_df <- df1 %>%
    full_join(df2, by = c("BIN1_CHR" = "BIN1_CHR", "BIN2_CHR" = "BIN2_CHR")) %>%
    filter(
      BIN1_START.x <= BIN1_END.y & BIN1_END.x >= BIN1_START.y,
      BIN2_START.x <= BIN2_END.y & BIN2_END.x >= BIN2_START.y
    )
  overlap_count <- nrow(overlap_df)
  
  
  overlap <- paste0("Overlap_", 1:overlap_count)
  unique_df1 <- paste0("DF1_Unique_", 1:(df1_count - overlap_count))
  unique_df2 <- paste0("DF2_Unique_", 1:(df2_count - overlap_count))
  vector1 <- c(overlap, unique_df1)
  vector2 <- c(overlap, unique_df2)
  
  venn_data <- list(vector1, vector2)
  names(venn_data) <- c(sample1, sample2)
  
  ggVennDiagram(venn_data, label="both", label_alpha = 0) +
    scale_fill_gradient(low="white", high="pink") +
    scale_color_manual(values = c("Arterioles"="skyblue2", "EDA"="salmon2")) +
    theme(legend.position = "None") + labs(title = plot_title) +
    coord_flip()
  
}


mc_NV_4kb <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/merged_NV_4kb_loops.tsv", header = T)
mc_NV_8kb <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/merged_NV_8kb_loops.tsv", header = T)
mc_NV_16kb <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/merged_NV_16kb_loops.tsv", header = T)
mc_VR_4kb <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/merged_VR_4kb_loops.tsv", header = T)
mc_VR_8kb <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/merged_VR_8kb_loops.tsv", header = T)
mc_VR_16kb <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/merged_VR_16kb_loops.tsv", header = T)

png(filename = "/xdisk/mliang1/qqiu/project/others/test/figure/fig2a.mc_4k.venn.png", 
    width = 220/96*300, height = 200/96*300, res = 300) 
venn_plot("mc_NV_4kb", "mc_VR_4kb")
dev.off()

png(filename = "/xdisk/mliang1/qqiu/project/others/test/figure/fig2a.mc_8k.venn.png", 
    width = 220/96*300, height = 200/96*300, res = 300) 
venn_plot("mc_NV_8kb", "mc_VR_8kb")
dev.off()

png(filename = "/xdisk/mliang1/qqiu/project/others/test/figure/fig2a.mc_16k.venn.png", 
    width = 220/96*300, height = 200/96*300, res = 300) 
venn_plot("mc_NV_16kb", "mc_VR_16kb")
dev.off()



cm_NV_10kb <- read.csv("/xdisk/mliang1/qqiu/project/others/test/data/NVCMC_newpara_10kb_intra.filtered.csv", header = T)
cm_NV_20kb <- read.csv("/xdisk/mliang1/qqiu/project/others/test/data/NVCMC_newpara_20kb_intra.filtered.csv", header = T)
cm_VR_10kb <- read.csv("/xdisk/mliang1/qqiu/project/others/test/data/VRCMC_newpara_10kb_intra.filtered.csv", header = T)
cm_VR_20kb <- read.csv("/xdisk/mliang1/qqiu/project/others/test/data/VRCMC_newpara_20kb_intra.filtered.csv", header = T)

png(filename = "/xdisk/mliang1/qqiu/project/others/test/figure/fig2b.cm_10k.venn.png", 
    width = 220/96*300, height = 200/96*300, res = 300) 
venn_plot("cm_NV_10kb", "cm_VR_10kb")
dev.off()

png(filename = "/xdisk/mliang1/qqiu/project/others/test/figure/fig2b.cm_20k.venn.png", 
    width = 220/96*300, height = 200/96*300, res = 300) 
venn_plot("cm_NV_20kb", "cm_VR_20kb")
dev.off()




################################################################################

regulatory_anno <- function(input_df, regulatory_ranges=regulatory_ranges){
  
  sample <- ifelse(grepl("NV", input_df), "Arterioles", "EDA")
  resolution <- strsplit(input_df, "_")[[1]][3] %>% gsub("kb", " kb", .)
  
  input_df <- standardize_column_names(get(input_df))
  
  mc_BIN1 <- GRanges(
    seqnames = input_df$BIN1_CHR,
    ranges = IRanges(start = input_df$BIN1_START, end = input_df$BIN1_END)
  )
  
  mc_BIN2 <- GRanges(
    seqnames = input_df$BIN2_CHR,
    ranges = IRanges(start = input_df$BIN2_START, end = input_df$BIN2_END)
  )
  
  overlap_BIN1 <- findOverlaps(mc_BIN1, regulatory_ranges)
  overlap_BIN2 <- findOverlaps(mc_BIN2, regulatory_ranges)
  
  input_df <- input_df %>%
    mutate(
      BIN1_overlap = countOverlaps(mc_BIN1, regulatory_ranges),
      BIN2_overlap = countOverlaps(mc_BIN2, regulatory_ranges),
      BIN1_regulatory_elements = sapply(seq_along(mc_BIN1), function(i) {
        if (any(queryHits(overlap_BIN1) == i)) {
          paste(unique(regulatory_ranges$V3[subjectHits(overlap_BIN1)[queryHits(overlap_BIN1) == i]]), collapse = ", ")
        } else {
          NA
        }
      }), 
      BIN2_regulatory_elements = sapply(seq_along(mc_BIN2), function(i) {
        if (any(queryHits(overlap_BIN2) == i)) {
          paste(unique(regulatory_ranges$V3[subjectHits(overlap_BIN2)[queryHits(overlap_BIN2) == i]]), collapse = ", ")
        } else {
          NA
        }
      }),
      regulatory_count = case_when(
        BIN1_overlap == 0 & BIN2_overlap == 0 ~ "None",
        BIN1_overlap != 0 & BIN2_overlap != 0 ~ "Both",
        .default = "One"
      ),
      sample = sample,
      resolution = resolution
    )
  return(input_df)
  
}


### apply all annotation for capture micro-c data
regulatory_ref <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/regulatory.reference.out", header = F, sep = "\t")
regulatory_ref$V1 <- paste0("chr", regulatory_ref$V1)

regulatory_ranges <- makeGRangesFromDataFrame(
  regulatory_ref,
  seqnames.field = "V1",
  start.field = "V4",
  end.field = "V5",
  keep.extra.columns = TRUE
)

mc_NV_4kb <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/merged_NV_4kb_loops.tsv", header = T)
mc_NV_8kb <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/merged_NV_8kb_loops.tsv", header = T)
mc_NV_16kb <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/merged_NV_16kb_loops.tsv", header = T)
mc_VR_4kb <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/merged_VR_4kb_loops.tsv", header = T)
mc_VR_8kb <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/merged_VR_8kb_loops.tsv", header = T)
mc_VR_16kb <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/merged_VR_16kb_loops.tsv", header = T)

mc_NV_4kb <- regulatory_anno("mc_NV_4kb", regulatory_ranges)
mc_NV_8kb <- regulatory_anno("mc_NV_8kb", regulatory_ranges)
mc_NV_16kb <- regulatory_anno("mc_NV_16kb", regulatory_ranges)
mc_VR_4kb <- regulatory_anno("mc_VR_4kb", regulatory_ranges)
mc_VR_8kb <- regulatory_anno("mc_VR_8kb", regulatory_ranges)
mc_VR_16kb <- regulatory_anno("mc_VR_16kb", regulatory_ranges)

plot_df <- rbind(mc_NV_4kb, mc_NV_8kb, mc_NV_16kb,
                mc_VR_4kb, mc_VR_8kb, mc_VR_16kb) %>%
  group_by(sample, resolution, regulatory_count) %>%
  summarize(Counts = n()) %>%
  mutate(resolution = factor(resolution, levels=c("4 kb", "8 kb", "16 kb")),
         regulatory_count = factor(regulatory_count, levels=c("None", "One", "Both")))

p <- ggplot(plot_df, aes(x=regulatory_count, y=Counts, fill=sample)) +
  geom_bar(stat="identity", width=.8, position = "dodge") +
  theme_classic() + 
  theme(axis.text = element_text(color="black", size=12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(color="black", size=12),
        strip.text = element_text(color="black", size=10)) +
  scale_fill_manual(values = c("Arterioles"="skyblue2", "EDA"="salmon2")) +
  labs(x = "Regulatory element occurrence across bins", fill = "Sample") +
  facet_wrap(~resolution)
print(p)
ggsave("/xdisk/mliang1/qqiu/project/others/test/figure/fig2c.reg_count.bar_plot.png", width=496/96, height=203/96, dpi=300)


both_prop <- c(
  sum(mc_NV_4kb$regulatory_count=="Both")/nrow(mc_NV_4kb),
  sum(mc_NV_8kb$regulatory_count=="Both")/nrow(mc_NV_8kb),
  sum(mc_NV_16kb$regulatory_count=="Both")/nrow(mc_NV_16kb),
  sum(mc_VR_4kb$regulatory_count=="Both")/nrow(mc_VR_4kb),
  sum(mc_VR_8kb$regulatory_count=="Both")/nrow(mc_VR_8kb),
  sum(mc_VR_16kb$regulatory_count=="Both")/nrow(mc_VR_16kb)
  )
mean(both_prop)
# [1] 0.91

one_prop <- c(
  sum(mc_NV_4kb$regulatory_count=="One")/nrow(mc_NV_4kb),
  sum(mc_NV_8kb$regulatory_count=="One")/nrow(mc_NV_8kb),
  sum(mc_NV_16kb$regulatory_count=="One")/nrow(mc_NV_16kb),
  sum(mc_VR_4kb$regulatory_count=="One")/nrow(mc_VR_4kb),
  sum(mc_VR_8kb$regulatory_count=="One")/nrow(mc_VR_8kb),
  sum(mc_VR_16kb$regulatory_count=="One")/nrow(mc_VR_16kb)
)
mean(one_prop)
# [1] 0.084



### apply non-promoter annotation for capture micro-c data
regulatory_ref <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/regulatory.reference.out", header = F, sep = "\t")
regulatory_ref <- regulatory_ref[regulatory_ref$V3!="promoter", c(1, 3:5)]
regulatory_ref$V1 <- paste0("chr", regulatory_ref$V1)

regulatory_ranges <- makeGRangesFromDataFrame(
  regulatory_ref,
  seqnames.field = "V1",
  start.field = "V4",
  end.field = "V5",
  keep.extra.columns = TRUE
)

cm_NV_10kb <- read.csv("/xdisk/mliang1/qqiu/project/others/test/data/NVCMC_newpara_10kb_intra.filtered.csv", header = T)
cm_NV_20kb <- read.csv("/xdisk/mliang1/qqiu/project/others/test/data/NVCMC_newpara_20kb_intra.filtered.csv", header = T)
cm_VR_10kb <- read.csv("/xdisk/mliang1/qqiu/project/others/test/data/VRCMC_newpara_10kb_intra.filtered.csv", header = T)
cm_VR_20kb <- read.csv("/xdisk/mliang1/qqiu/project/others/test/data/VRCMC_newpara_20kb_intra.filtered.csv", header = T)

cm_NV_10kb <- regulatory_anno("cm_NV_10kb", regulatory_ranges)
cm_NV_20kb <- regulatory_anno("cm_NV_20kb", regulatory_ranges)
cm_VR_10kb <- regulatory_anno("cm_VR_10kb", regulatory_ranges)
cm_VR_20kb <- regulatory_anno("cm_VR_20kb", regulatory_ranges)

plot_df <- rbind(cm_NV_10kb, cm_NV_20kb,
                cm_VR_10kb, cm_VR_20kb) %>%
  group_by(sample, resolution, regulatory_count) %>%
  summarize(Counts = n()) %>%
  mutate(resolution = factor(resolution, levels=c("10 kb", "20 kb")),
         regulatory_count = factor(regulatory_count, levels=c("None", "One", "Both")))
p <- ggplot(plot_df, aes(x=regulatory_count, y=Counts, fill=sample)) +
  geom_bar(stat="identity", width=.8, position = "dodge") +
  theme_classic() + 
  theme(axis.text = element_text(color="black", size=12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(color="black", size=12),
        strip.text = element_text(color="black", size=10)) +
  scale_fill_manual(values = c("Arterioles"="skyblue2", "EDA"="salmon2")) +
  labs(x = "Regulatory element occurrence across bins", fill = "Sample") +
  facet_wrap(~resolution)
print(p)
ggsave("/xdisk/mliang1/qqiu/project/others/test/figure/figs8.reg_count.bar_plot.png", width=496/96, height=203/96, dpi=300)



################################################################################

regulatory_anno_expd <- function(input_df, regulatory_ranges=regulatory_ranges){
  
  sample <- ifelse(grepl("NV", input_df), "Arterioles", "EDA")
  resolution <- strsplit(input_df, "_")[[1]][3] %>% gsub("kb", " kb", .)
  
  input_df <- standardize_column_names(get(input_df))
  
  mc_BIN1 <- GRanges(
    seqnames = input_df$BIN1_CHR,
    ranges = IRanges(start = input_df$BIN1_START, end = input_df$BIN1_END)
  )
  
  mc_BIN2 <- GRanges(
    seqnames = input_df$BIN2_CHR,
    ranges = IRanges(start = input_df$BIN2_START, end = input_df$BIN2_END)
  )
  
  overlap_BIN1 <- findOverlaps(mc_BIN1, regulatory_ranges)
  overlap_BIN2 <- findOverlaps(mc_BIN2, regulatory_ranges)
  
  expanded_BIN1 <- data.frame(
    BIN_ID = queryHits(overlap_BIN1),
    BIN1_regulatory_element = regulatory_ranges$regulatory_element[subjectHits(overlap_BIN1)],
    BIN1_transcript = regulatory_ranges$transcript_id[subjectHits(overlap_BIN1)],
    BIN1_gene = regulatory_ranges$gene_id[subjectHits(overlap_BIN1)]
  ) %>% 
    mutate(
      BIN1_reg_abb = toupper(substr(BIN1_regulatory_element, 1, 1)),
      BIN1_reg_abb = ifelse(is.na(BIN1_reg_abb), "", BIN1_reg_abb)
    )
  
  expanded_BIN2 <- data.frame(
    BIN_ID = queryHits(overlap_BIN2),
    BIN2_regulatory_element = regulatory_ranges$regulatory_element[subjectHits(overlap_BIN2)],
    BIN2_transcript = regulatory_ranges$transcript_id[subjectHits(overlap_BIN2)],
    BIN2_gene = regulatory_ranges$gene_id[subjectHits(overlap_BIN2)]
  ) %>%
    mutate(
      BIN2_reg_abb = toupper(substr(BIN2_regulatory_element, 1, 1)),
      BIN2_reg_abb = ifelse(is.na(BIN2_reg_abb), "", BIN2_reg_abb)
    )
  
  expanded_input_df <- input_df %>%
    mutate(BIN_ID = row_number()) %>%
    left_join(expanded_BIN1, by = c("BIN_ID")) %>%
    left_join(expanded_BIN2, by = c("BIN_ID")) %>%
    mutate(
      sample = sample,
      resolution = resolution,
      interaction_cat = sapply(1:nrow(.), function(i) {
        paste0(sort(c(BIN1_reg_abb[i], BIN2_reg_abb[i])), collapse = "")
      })
    ) %>% unique()
    
  return(expanded_input_df)
  
}


###
canonical_chromosomes <- c(paste0("chr", 1:22), "chrX", "chrY")

ref <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/Homo_sapiens.GRCh38.112.transcript.gtf", header = F, sep = '\t')
ref$transcript_id <- str_extract(ref$V9, "transcript_id [^;]+") %>%
  gsub("transcript_id ", "", .)
ref$gene_id <- str_extract(ref$V9, "gene_id [^;]+") %>%
  gsub("gene_id ", "", .)

transcript_ranges <- GRanges(
  seqnames = Rle(paste0("chr", ref$V1)),
  ranges = IRanges(start = ref$V4, end = ref$V5),
  strand = Rle(ref$V7),
  transcript_id = ref$transcript_id,
  gene_id = ref$gene_id
)

regulatory_ref = read.table("/xdisk/mliang1/qqiu/project/others/test/data/regulatory.reference.out", header = F, sep = "\t")
regulatory_ref = regulatory_ref[regulatory_ref$V3 %in% c("TF_binding_site", "enhancer", "promoter"), c(1, 3:5)]
regulatory_ref$V1 = paste0("chr", regulatory_ref$V1)
regulatory_ref$reg_idx <- seq_len(nrow(regulatory_ref))

regulatory_ranges <- makeGRangesFromDataFrame(
  regulatory_ref,
  seqnames.field = "V1",
  start.field = "V4",
  end.field = "V5",
  keep.extra.columns = TRUE
)
regulatory_ranges <- regulatory_ranges[seqnames(regulatory_ranges) %in% canonical_chromosomes]
regulatory_ranges <- keepSeqlevels(regulatory_ranges, canonical_chromosomes, pruning.mode = "coarse")

extended_regulatory_ranges <- resize(regulatory_ranges, 
                                     width(regulatory_ranges) + 10000,  # extend by 5kb on each side
                                     fix = "center")

overlaps <- findOverlaps(extended_regulatory_ranges, transcript_ranges)
overlap_df <- data.frame(
  reg_idx = extended_regulatory_ranges$reg_idx[queryHits(overlaps)],
  transcript_id = transcript_ranges$transcript_id[subjectHits(overlaps)],
  gene_id = transcript_ranges$gene_id[subjectHits(overlaps)]
)

all_regulatory_df <- data.frame(
  seqnames = seqnames(regulatory_ranges),
  start = start(regulatory_ranges),
  end = end(regulatory_ranges),
  regulatory_element = regulatory_ranges$V3,
  reg_idx = regulatory_ranges$reg_idx
)

final_df <- merge(all_regulatory_df, overlap_df, 
                  by = c("reg_idx"), 
                  all.x = TRUE)
final_df[is.na(final_df)] <- ""

regulatory_ranges <- GRanges(
  seqnames = final_df$seqnames,
  ranges = IRanges(start = final_df$start, 
                   end = final_df$end),
  regulatory_id = paste(final_df$seqnames, final_df$start, final_df$end, sep = "_"),
  regulatory_element = final_df$regulatory_element,
  transcript_id = final_df$transcript_id,
  gene_id = final_df$gene_id
)


mc_NV_4kb <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/merged_NV_4kb_loops.tsv", header = T)
mc_NV_8kb <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/merged_NV_8kb_loops.tsv", header = T)
mc_NV_16kb <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/merged_NV_16kb_loops.tsv", header = T)
mc_VR_4kb <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/merged_VR_4kb_loops.tsv", header = T)
mc_VR_8kb <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/merged_VR_8kb_loops.tsv", header = T)
mc_VR_16kb <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/merged_VR_16kb_loops.tsv", header = T)

mc_NV_4kb <- regulatory_anno_expd("mc_NV_4kb", regulatory_ranges)
mc_NV_8kb <- regulatory_anno_expd("mc_NV_8kb", regulatory_ranges)
mc_NV_16kb <- regulatory_anno_expd("mc_NV_16kb", regulatory_ranges)
mc_VR_4kb <- regulatory_anno_expd("mc_VR_4kb", regulatory_ranges)
mc_VR_8kb <- regulatory_anno_expd("mc_VR_8kb", regulatory_ranges)
mc_VR_16kb <- regulatory_anno_expd("mc_VR_16kb", regulatory_ranges)


plot_df <- rbind(mc_NV_4kb, mc_NV_8kb, mc_NV_16kb,
                mc_VR_4kb, mc_VR_8kb, mc_VR_16kb) %>%
  dplyr::select(sample, resolution, BIN_ID, interaction_cat) %>% unique() %>%
  group_by(sample, resolution, interaction_cat) %>%
  summarize(Counts = n()) %>%
  mutate(resolution = factor(resolution, levels=c("4 kb", "8 kb", "16 kb")),
         interaction_cat = factor(interaction_cat, levels=c("EE", "EP", "ET", "PP", "PT", "TT"))) %>%
  filter(!is.na(interaction_cat))

p <- ggplot(plot_df, aes(y=interaction_cat, x=Counts, fill=sample)) +
  geom_bar(stat="identity", width=.8, position = "dodge") +
  theme_classic() + 
  theme(axis.text = element_text(color="black", size=12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(color="black", size=12),
        strip.text = element_text(color="black", size=10)) +
  scale_fill_manual(values = c("Arterioles"="skyblue2", "EDA"="salmon2")) +
  labs(x = "Counts", y = "Interaction of regulatory elements", fill = "Sample") +
  facet_wrap(~resolution, scales = "free_x")
print(p)
ggsave("/xdisk/mliang1/qqiu/project/others/test/figure/figs7.mc.int_reg_count.bar_plot.png", width=496/96, height=340/96, dpi=300)




### build updated regulatory ranges for capture micro-c data
canonical_chromosomes <- c(paste0("chr", 1:22), "chrX", "chrY")

ref <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/Homo_sapiens.GRCh38.112.transcript.gtf", header = F, sep = '\t')
ref$transcript_id <- str_extract(ref$V9, "transcript_id [^;]+") %>%
  gsub("transcript_id ", "", .)
ref$gene_id <- str_extract(ref$V9, "gene_id [^;]+") %>%
  gsub("gene_id ", "", .)

transcript_ranges <- GRanges(
  seqnames = Rle(paste0("chr", ref$V1)),
  ranges = IRanges(start = ref$V4, end = ref$V5),
  strand = Rle(ref$V7),
  transcript_id = ref$transcript_id,
  gene_id = ref$gene_id
)
promoter_ranges <- promoters(transcript_ranges, upstream = 1000, downstream = 500)
mcols(promoter_ranges)$regulatory_id <- paste(seqnames(promoter_ranges), start(promoter_ranges), end(promoter_ranges), sep = "_")
mcols(promoter_ranges)$regulatory_element <- "promoter"
mcols(promoter_ranges)$transcript_id <- transcript_ranges$transcript_id
mcols(promoter_ranges)$gene_id <- transcript_ranges$gene_id
promoter_ranges <- promoter_ranges[seqnames(promoter_ranges) %in% canonical_chromosomes]
promoter_ranges <- keepSeqlevels(promoter_ranges, canonical_chromosomes, pruning.mode = "coarse")


regulatory_ref <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/regulatory.reference.out", header = F, sep = "\t")
regulatory_ref <- regulatory_ref[regulatory_ref$V3 %in% c("TF_binding_site", "enhancer"), c(1, 3:5)]
regulatory_ref$V1 <- paste0("chr", regulatory_ref$V1)
regulatory_ref$reg_idx <- seq_len(nrow(regulatory_ref))

regulatory_ranges <- makeGRangesFromDataFrame(
  regulatory_ref,
  seqnames.field = "V1",
  start.field = "V4",
  end.field = "V5",
  keep.extra.columns = TRUE
)
regulatory_ranges <- regulatory_ranges[seqnames(regulatory_ranges) %in% canonical_chromosomes]
regulatory_ranges <- keepSeqlevels(regulatory_ranges, canonical_chromosomes, pruning.mode = "coarse")

extended_regulatory_ranges <- resize(regulatory_ranges, 
                                     width(regulatory_ranges) + 10000,  # extend by 5kb on each side
                                     fix = "center")

overlaps <- findOverlaps(extended_regulatory_ranges, transcript_ranges)
overlap_df <- data.frame(
  reg_idx = extended_regulatory_ranges$reg_idx[queryHits(overlaps)],
  transcript_id = transcript_ranges$transcript_id[subjectHits(overlaps)],
  gene_id = transcript_ranges$gene_id[subjectHits(overlaps)]
)

all_regulatory_df <- data.frame(
  seqnames = seqnames(regulatory_ranges),
  start = start(regulatory_ranges),
  end = end(regulatory_ranges),
  regulatory_element = regulatory_ranges$V3,
  reg_idx = regulatory_ranges$reg_idx
)

final_df <- merge(all_regulatory_df, overlap_df, 
                  by = c("reg_idx"), 
                  all.x = TRUE)
final_df[is.na(final_df)] <- ""

final_granges <- GRanges(
  seqnames = final_df$seqnames,
  ranges = IRanges(start = final_df$start, 
                   end = final_df$end),
  regulatory_id = paste(final_df$seqnames, final_df$start, final_df$end, sep = "_"),
  regulatory_element = final_df$regulatory_element,
  transcript_id = final_df$transcript_id,
  gene_id = final_df$gene_id
)

regulatory_ranges = c(final_granges, promoter_ranges)


cm_NV_10kb <- read.csv("/xdisk/mliang1/qqiu/project/others/test/data/NVCMC_newpara_10kb_intra.filtered.csv", header = T)
cm_NV_20kb <- read.csv("/xdisk/mliang1/qqiu/project/others/test/data/NVCMC_newpara_20kb_intra.filtered.csv", header = T)
cm_VR_10kb <- read.csv("/xdisk/mliang1/qqiu/project/others/test/data/VRCMC_newpara_10kb_intra.filtered.csv", header = T)
cm_VR_20kb <- read.csv("/xdisk/mliang1/qqiu/project/others/test/data/VRCMC_newpara_20kb_intra.filtered.csv", header = T)

cm_NV_10kb <- regulatory_anno_expd("cm_NV_10kb", regulatory_ranges)
cm_NV_20kb <- regulatory_anno_expd("cm_NV_20kb", regulatory_ranges)
cm_VR_10kb <- regulatory_anno_expd("cm_VR_10kb", regulatory_ranges)
cm_VR_20kb <- regulatory_anno_expd("cm_VR_20kb", regulatory_ranges)


plot_df <- rbind(cm_NV_10kb, cm_NV_20kb, 
                cm_VR_10kb, cm_VR_20kb) %>%
  dplyr::select(sample, resolution, BIN_ID, interaction_cat) %>% unique() %>%
  group_by(sample, resolution, interaction_cat) %>%
  summarize(Counts = n()) %>%
  mutate(resolution = factor(resolution, levels=c("10 kb", "20 kb")),
         interaction_cat = factor(interaction_cat, levels=c("EE", "EP", "ET", "PP", "PT", "TT"))) %>%
  filter(!is.na(interaction_cat))

p <- ggplot(plot_df, aes(y=interaction_cat, x=Counts, fill=sample)) +
  geom_bar(stat="identity", width=.8, position = "dodge") +
  theme_classic() + 
  theme(axis.text = element_text(color="black", size=12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(color="black", size=12),
        strip.text = element_text(color="black", size=10)) +
  scale_fill_manual(values = c("Arterioles"="skyblue2", "EDA"="salmon2")) +
  labs(x = "Counts", y = "Interaction of regulatory elements", fill = "Sample") +
  facet_wrap(~resolution, scales = "free_x")
print(p)
ggsave("/xdisk/mliang1/qqiu/project/others/test/figure/figs7.cm.int_reg_count.bar_plot.png", width=496/96, height=340/96, dpi=300)



save(mc_NV_4kb, mc_NV_8kb, mc_NV_16kb,
     mc_VR_4kb, mc_VR_8kb, mc_VR_16kb, 
     cm_NV_10kb, cm_NV_20kb, 
     cm_VR_10kb, cm_VR_20kb, file = "/xdisk/mliang1/qqiu/project/others/test/data/mc_cm.regulatory_anno.RData")








gene_expd <- function(input_df, filtered_regulatory_ranges=filtered_regulatory_ranges, expr=expr, meth=meth){
  
  sample <- ifelse(grepl("NV", input_df), "Arterioles", "EDA")
  resolution <- strsplit(input_df, "_")[[1]][3] %>% gsub("kb", " kb", .)
  
  if(sample=="Arterioles"){
    expr_list <- colnames(expr)[grep("namR$", colnames(expr))]
    meth_list <- colnames(meth)[grep("^HNV", colnames(meth))]
  }else{
    expr_list <- colnames(expr)[grep("vrmR$", colnames(expr))]
    meth_list <- colnames(meth)[grep("^HVR", colnames(meth))]
  }
  expr_avg <- data.frame(gene_id = rownames(expr),
                         expr = rowMeans(log2(expr[, expr_list]+1)))
  meth_avg <- data.frame(gene_id = rownames(meth),
                         meth = rowMeans(meth[, meth_list]))
  
  input_df <- standardize_column_names(get(input_df))
  input_df <- unique(rbind(
    setNames(input_df[, c("BIN1_gene", "interaction_cat")], c("gene_id", "interaction_cat")),
    setNames(input_df[, c("BIN2_gene", "interaction_cat")], c("gene_id", "interaction_cat"))
  )) %>% filter(!is.na(gene_id)) %>% mutate(interaction_type="")
  
  gene_df = unique(mcols(filtered_regulatory_ranges)[, c("regulatory_element", "gene_id")]) %>%
    as.data.frame() %>%
    filter(!(gene_id %in% input_df$gene_id)) %>%
    mutate(interaction_cat = toupper(substr(regulatory_element, 1, 1)),
           interaction_type = "Control.") %>%
    dplyr::select(gene_id, interaction_cat, interaction_type) %>%
    unique()
  
  merged_df = rbind(input_df, gene_df) %>%
    mutate(interaction_type = paste0(interaction_type, interaction_cat)) %>%
    left_join(expr_avg, by = c("gene_id")) %>%
    left_join(meth_avg, by = c("gene_id")) %>%
    mutate(
      sample = sample,
      resolution = resolution
    )
  
  return(merged_df)
  
}

filtered_regulatory_ranges <- regulatory_ranges

expr = read.csv("/xdisk/mliang1/qqiu/project/others/test/data/vessel.gene.FPKM.rm_013.csv", header = T, row.names = 1)
meth = read.csv("/xdisk/mliang1/qqiu/project/others/test/data/vessel.RRBS.dat.cln.promoter.methyl.out", header = T, row.names = 1)

load("/xdisk/mliang1/qqiu/project/others/test/data/mc_cm.regulatory_anno.RData")

mc_NV_4kb_trans <- gene_expd("mc_NV_4kb", filtered_regulatory_ranges, expr, meth)
mc_NV_8kb_trans <- gene_expd("mc_NV_8kb", filtered_regulatory_ranges, expr, meth)
mc_NV_16kb_trans <- gene_expd("mc_NV_16kb", filtered_regulatory_ranges, expr, meth)
mc_VR_4kb_trans <- gene_expd("mc_VR_4kb", filtered_regulatory_ranges, expr, meth)
mc_VR_8kb_trans <- gene_expd("mc_VR_8kb", filtered_regulatory_ranges, expr, meth)
mc_VR_16kb_trans <- gene_expd("mc_VR_16kb", filtered_regulatory_ranges, expr, meth)


plot_df <- rbind(mc_NV_4kb_trans, mc_NV_8kb_trans, mc_NV_16kb_trans,
                mc_VR_4kb_trans, mc_VR_8kb_trans, mc_VR_16kb_trans) %>%
  mutate(resolution = factor(resolution, levels=c("4 kb", "8 kb", "16 kb")),
         interaction_type = factor(interaction_type)) %>%
  filter(!is.na(interaction_type) & !(interaction_type %in% c("E", "P", "T")))

median_values <- plot_df %>%
  filter(grepl("Control", interaction_type)) %>%
  group_by(resolution) %>%
  summarize(median_expr = median(expr, na.rm = TRUE))

p <- ggplot(plot_df, aes(y=interaction_type, x=expr, fill=sample)) +
  geom_boxplot(width=.8, position = "dodge", outlier.shape = NA) +
  theme_classic() + 
  theme(axis.text = element_text(color="black", size=12),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(color="black", size=12),
        strip.text = element_text(color="black", size=10)) +
  scale_fill_manual(values = c("Arterioles"="skyblue2", "EDA"="salmon2")) +
  labs(x = "Gene expression (log2(FPKM+1))", y = "Interaction of regulatory elements", fill = "Sample") +
  geom_vline(data = median_values, aes(xintercept = median_expr), color = "red", linetype = "dashed") +
  facet_wrap(~resolution)
print(p)
ggsave("/xdisk/mliang1/qqiu/project/others/test/figure/fig2e.mc.int_reg_expr.boxplot.png", width=496/96, height=340/96, dpi=300)



probe_trans <- read.csv("/xdisk/mliang1/qqiu/project/others/test/data/h_probes_v1.0.trans_gene.out", header = T, sep = ',')
filtered_regulatory_ranges <- regulatory_ranges[regulatory_ranges$transcript_id %in% probe_trans$transcript_id, ]

cm_NV_10kb_trans <- gene_expd("cm_NV_10kb", filtered_regulatory_ranges, expr, meth)
cm_NV_20kb_trans <- gene_expd("cm_NV_20kb", filtered_regulatory_ranges, expr, meth)
cm_VR_10kb_trans <- gene_expd("cm_VR_10kb", filtered_regulatory_ranges, expr, meth)
cm_VR_20kb_trans <- gene_expd("cm_VR_20kb", filtered_regulatory_ranges, expr, meth)

plot_df <- rbind(cm_NV_10kb_trans, cm_NV_20kb_trans,
                cm_VR_10kb_trans, cm_VR_20kb_trans) %>%
  mutate(resolution = factor(resolution, levels=c("10 kb", "20 kb")),
         interaction_type = factor(interaction_type)) %>%
  filter(!is.na(interaction_type) & !(interaction_type %in% c("E", "P", "T")))

p <- ggplot(plot_df, aes(y=interaction_type, x=expr, fill=sample)) +
  geom_boxplot(width=.8, position = "dodge", outlier.shape = NA) +
  theme_classic() + 
  theme(axis.text = element_text(color="black", size=12),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(color="black", size=12),
        strip.text = element_text(color="black", size=10)) +
  scale_fill_manual(values = c("Arterioles"="skyblue2", "EDA"="salmon2")) +
  labs(x = "Gene expression (log2(FPKM+1))", y = "Interaction of regulatory elements", fill = "Sample") +
  facet_wrap(~resolution)
print(p)
ggsave("/xdisk/mliang1/qqiu/project/others/test/figure/fig2e.cm.int_reg_expr.boxplot.png", width=496/96, height=340/96, dpi=300)






################################################################################

promoter_count <- function(input_df, transcript_ranges=transcript_ranges, expr=expr){
  
  sample <- ifelse(grepl("NV", input_df), "Arterioles", "EDA")
  resolution <- strsplit(input_df, "_")[[1]][3] %>% gsub("kb", " kb", .)
  
  if(sample=="Arterioles"){
    expr_list <- colnames(expr)[grep("namR$", colnames(expr))]
    meth_list <- colnames(meth)[grep("^HNV", colnames(meth))]
  }else{
    expr_list <- colnames(expr)[grep("vrmR$", colnames(expr))]
    meth_list <- colnames(meth)[grep("^HVR", colnames(meth))]
  }
  expr_avg <- data.frame(gene_id = rownames(expr),
                         expr = rowMeans(log2(expr[, expr_list]+1)))
  meth_avg <- data.frame(gene_id = rownames(meth),
                         meth = rowMeans(meth[, meth_list]))
  
  input_df <- standardize_column_names(get(input_df))
  
  mc_BIN1 <- GRanges(
    seqnames = input_df$BIN1_CHR,
    ranges = IRanges(start = input_df$BIN1_START, end = input_df$BIN1_END)
  )
  
  mc_BIN2 <- GRanges(
    seqnames = input_df$BIN2_CHR,
    ranges = IRanges(start = input_df$BIN2_START, end = input_df$BIN2_END)
  )
  
  overlap_BIN1 <- findOverlaps(mc_BIN1, transcript_ranges)
  overlap_BIN2 <- findOverlaps(mc_BIN2, transcript_ranges)
  
  trans_BIN1 <- data.frame(
    BIN_id = queryHits(overlap_BIN1),
    gene_id = transcript_ranges$gene_id[subjectHits(overlap_BIN1)]
  ) %>% unique() %>% group_by(gene_id) %>% 
    summarise(count = n())
    
  trans_BIN2 <- data.frame(
    BIN_id = queryHits(overlap_BIN2),
    gene_id = transcript_ranges$gene_id[subjectHits(overlap_BIN2)]
  ) %>% unique() %>% group_by(gene_id) %>% 
    summarise(count = n())
  
  trans_expd = data.frame(gene_id = unique(transcript_ranges$gene_id)) %>%
    left_join(trans_BIN1, by = c("gene_id")) %>%
    left_join(trans_BIN2, by = c("gene_id")) %>%
    mutate(across(everything(), ~replace_na(., 0)),
           count = count.x + count.y) %>%
    left_join(expr_avg, by = c("gene_id")) %>%
    left_join(meth_avg, by = c("gene_id")) %>%
    mutate(
      sample = sample,
      resolution = resolution
    ) %>% dplyr::select(-count.x, -count.y) %>% unique()
  
  return(trans_expd)
  
}

ref <- read.table("/xdisk/mliang1/qqiu/project/others/test/data/Homo_sapiens.GRCh38.112.transcript.gtf", header = F, sep = '\t')
probe_trans <- read.csv("/xdisk/mliang1/qqiu/project/others/test/data/h_probes_v1.0.trans_gene.out", header = T, sep = ',')
ref$transcript_id <- str_extract(ref$V9, "transcript_id [^;]+") %>%
  gsub("transcript_id ", "", .)
ref$gene_id <- str_extract(ref$V9, "gene_id [^;]+") %>%
  gsub("gene_id ", "", .)
ref = ref[ref$transcript_id %in% probe_trans$transcript_id, ]

transcript_ranges <- GRanges(
  seqnames = Rle(paste0("chr", ref$V1)),
  ranges = IRanges(start = ref$V4, end = ref$V5),
  strand = Rle(ref$V7),
  transcript_id = ref$transcript_id,
  gene_id = ref$gene_id
)

expr = read.csv("/xdisk/mliang1/qqiu/project/others/test/data/vessel.gene.FPKM.rm_013.csv", header = T, row.names = 1)
meth = read.csv("/xdisk/mliang1/qqiu/project/others/test/data/vessel.RRBS.dat.cln.promoter.methyl.out", header = T, row.names = 1)


cm_NV_10kb <- read.csv("/xdisk/mliang1/qqiu/project/others/test/data/NVCMC_newpara_10kb_intra.filtered.csv", header = T)
cm_NV_20kb <- read.csv("/xdisk/mliang1/qqiu/project/others/test/data/NVCMC_newpara_20kb_intra.filtered.csv", header = T)
cm_VR_10kb <- read.csv("/xdisk/mliang1/qqiu/project/others/test/data/VRCMC_newpara_10kb_intra.filtered.csv", header = T)
cm_VR_20kb <- read.csv("/xdisk/mliang1/qqiu/project/others/test/data/VRCMC_newpara_20kb_intra.filtered.csv", header = T)

cm_NV_10kb_trans <- promoter_count("cm_NV_10kb", transcript_ranges, expr)
cm_NV_20kb_trans <- promoter_count("cm_NV_20kb", transcript_ranges, expr)
cm_VR_10kb_trans <- promoter_count("cm_VR_10kb", transcript_ranges, expr)
cm_VR_20kb_trans <- promoter_count("cm_VR_20kb", transcript_ranges, expr)


boxplot_mod <- function(input_df){
  
  sample <- ifelse(grepl("NV", input_df), "Arterioles", "EDA")
  resolution <- strsplit(input_df, "_")[[1]][3] %>% gsub("kb", " kb", .)
  plot_title <- paste0(sample, " (", resolution, ")\n")
  
  input_df = get(input_df) %>%
    mutate(resolution = factor(resolution, levels=c("10 kb", "20 kb")),
           count_fct = ifelse(count > 50, ">50", as.character(count)),
           count_fct = factor(count_fct, levels = c(as.character(0:50), ">50")))
  
  lm = coef(lm(expr ~ count, data = input_df)); intercept=as.numeric(lm[1]); slope=as.numeric(lm[2])
  cor=cor.test(input_df$expr,as.numeric(input_df$count), method = "spearman", exact = FALSE); cor_est=cor$estimate; cor_p=cor$p.value 
  cor_label = paste0("R = ", round(cor_est, 2), ", p-value = ", formatC(cor_p, format = "e", digits = 2))
  p = ggplot(input_df, aes(y=expr, x=count_fct, fill=sample)) +
    geom_boxplot(outlier.shape=NA) +
    theme_classic() + 
    theme(axis.text = element_text(color="black", size=10),
          axis.title = element_text(color="black", size=12),
          legend.position = "none") +
    scale_fill_manual(values = c("Arterioles"="skyblue2", "EDA"="salmon2")) +
    labs(x = "Promoter interaction count", y = "Gene expression (log2(FPKM+1))", title = plot_title) +
    scale_x_discrete(labels = function(x) {
      ifelse(x %in% c("0", "10", "20", "30", "40", ">50"), x, "")
    }) +
    annotate(geom="text", x=max(as.numeric(input_df$count_fct), na.rm = TRUE)*0.6, y=max(input_df$expr, na.rm = TRUE)*0.95, 
             label = cor_label, color = "red", size = 4.5) +
    geom_abline(intercept = intercept, slope = slope, color = "red", linetype = "dashed")
  
  print(p)

}

boxplot_mod("cm_NV_20kb_trans")
ggsave("/xdisk/mliang1/qqiu/project/others/test/figure/fig2f.cm_NV_20kb_promoter_count_expr_cor.png", width=381/96, height=327/96, dpi=300)

boxplot_mod("cm_VR_20kb_trans")
ggsave("/xdisk/mliang1/qqiu/project/others/test/figure/fig2f.cm_VR_20kb_promoter_count_expr_cor.png", width=381/96, height=327/96, dpi=300)




