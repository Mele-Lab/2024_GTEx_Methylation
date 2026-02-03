
# plot CV values for different K values in admixture

library(tidyverse)


first_dir <- "/Users/mariasopenar/cluster/"

cv <- readLines(paste0(first_dir, "Projects/GTEx_v8/Methylation/cv_errors.txt")) %>%
  tibble(line = .) %>%
  extract(
    line,
    into = c("K", "CV"),
    regex = "CV error \\(K=([0-9]+)\\): ([0-9.]+)"
  ) %>%
  mutate(
    K = as.integer(K),
    CV = as.numeric(CV)
  ) %>%
  arrange(K)

pdf(paste0(first_dir, "Projects/GTEx_v8/Methylation/Plots/","Admixture_cross_validation.pdf"), width = 3.74, height = 3.74)
ggplot(cv, aes(x = as.factor(K), y = CV)) +
  geom_point(size = 4, color="#F0AE21") +
  theme_classic() +
  labs(
    x = "Number of ancestral populations (K)",
    y = "Cross-validation error",
  )+ theme(legend.title = element_blank(),
           axis.text.x = element_text(colour="black", size=11),
           axis.text.y = element_text(colour="black", size=11),
           legend.text = element_text(colour="black", size=13),
           axis.title = element_text(size=13),
           legend.spacing.y = unit(-0.05, "cm"),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.border = element_rect(colour = "black", linewidth=1))
dev.off()

