


md_info <- read.delim(paste0(first_dir, "Projects/GTEx_v8/Methylation/Data/GTEx_Sample_Attributes.GRU.txt"),
                           header = TRUE,
                           row.names = 1,
                           quote = "",
                           fill = TRUE)


x <- read.table(paste0(first_dir, "Projects/GTEx_v8/Methylation/GTEx_Analysis_v8_xCell_scores_7_celltypes.txt"), header=T)
metadata <- readRDS(paste0(project_path, 'Tissues/', tissue, "/metadata.rds"))


library(dplyr)
library(stringr)

md_info_short <- md_info %>%
  mutate(
    SUBJID = str_extract(SAMPID, "^[^-]+-[^-]+"),
    TISSUE = SMTSD %>%
      str_remove("\\s*-.*$") %>%
      str_trim() %>%
      str_replace_all("\\s+", ""),   
      SUBJID_XCell = str_replace_all(SAMPID, "-", ".")
  ) %>% select( SUBJID, SAMPID, TISSUE, SUBJID_XCell)

saveRDS(md_info_short, paste0(first_dir, "Projects/GTEx_v8/Methylation/Data/donor_sample_correspondance.rds"))
setdiff(colnames(x), md_info_short$SUBJID_XCell)



cell_props_all <- read.table(paste0(first_dir, "Projects/GTEx_v8/Methylation/GTEx_Analysis_v8_xCell_scores_7_celltypes.txt"), header=T, row.names = 1)
md_info_short <- readRDS(paste0(first_dir, "Projects/GTEx_v8/Methylation/Data/donor_sample_correspondance.rds"))

md_info_short$TISSUE <- gsub("Colon", "ColonTransverse", md_info_short$TISSUE)
md_info_short$TISSUE <- gsub("Breast", "BreastMammaryTissue", md_info_short$TISSUE)
md_info_short$TISSUE <- gsub("Muscle", "MuscleSkeletal", md_info_short$TISSUE)
md_info_short$TISSUE <- gsub("Kidney", "KidneyCortex", md_info_short$TISSUE)
tissues <- c("BreastMammaryTissue", "ColonTransverse" ,"KidneyCortex", "Lung", "MuscleSkeletal" ,"Ovary", "Prostate", "Testis", "WholeBlood")




get_stats_XCell <- function(tissue){
  metadata <- readRDS(paste0(project_path, 'Tissues/', tissue, "/metadata.rds"))
  donor_sample <- md_info_short[md_info_short$TISSUE == tissue,]
  cell_props <- t(cell_props_all[, colnames(cell_props_all) %in% unique(donor_sample$SUBJID_XCell)])
  rownames(cell_props) <- str_extract(rownames(cell_props), "^[^.]+\\.[^.]+")
  rownames(cell_props) <- gsub("\\.", "-",  rownames(cell_props)) 
  length(metadata$SUBJID)
  length(intersect(rownames(cell_props), metadata$SUBJID))
  print("Correcting per cell type composition based on gene expression ---")
  print(paste0("We keep ",  length(intersect(rownames(cell_props), metadata$SUBJID)), " samples out of ",     length(metadata$SUBJID)))
  stats <- data.frame(c(length(intersect(rownames(cell_props), metadata$SUBJID))), length(metadata$SUBJID))
  rownames(stats) <- tissue
  colnames(stats)<- c("N_donors_XCell", "N_donors_total")
  return(stats)
}

stats <- lapply(tissues, function(tissue) get_stats_XCell(tissue))
stats_df <- do.call(rbind.data.frame,stats)


library(tidyr)
plot_df <- stats_df %>%
  tibble::rownames_to_column("TISSUE") %>%
  mutate(
    N_not_XCell = pmax(N_donors_total - N_donors_XCell, 0)  # safety
  ) %>%
  select(TISSUE, XCell_Scores = N_donors_XCell, Not_XCell_Scores = N_not_XCell) %>%
  pivot_longer(cols = c(XCell_Scores, Not_XCell_Scores), names_to = "Status", values_to = "N") %>%
  # order tissues by total donors (optional)
  group_by(TISSUE) %>%
  mutate(Total = sum(N)) %>%
  ungroup() %>%
  mutate(TISSUE = reorder(TISSUE, Total))

ggplot(plot_df, aes(x = TISSUE, y = N, alpha = Status)) +
  geom_col() +
  coord_flip() +
  labs(x = NULL, y = "Number of donors", fill = NULL) +
  theme_bw()+xlab("Tissue") +scale_alpha_manual(values = c(0.2, 1), name="")+
  theme(#legend.position = "none",
    axis.text.x = element_text(colour="black", size=11),
    axis.text.y = element_text(colour="black", size=13),
    legend.text = element_text(colour="black", size=12),
    axis.title.x = element_text(size=13),
    legend.spacing.y = unit(-0.05, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", linewidth=1)) 


