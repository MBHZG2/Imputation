#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
db <- read.table(args[1])
colnames(db) <- c("chr", "SNP", "POS", "A1", "A2", "Q", "PASS", "DR2", "MAF")

# Install and load the required packages if not already installed
required_packages <- c("dplyr", "ggplot2")
for (package in required_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
}

# Load the required packages
library(ggplot2)
library(dplyr)

# Convert DR2 to numeric
db$DR2 <- as.numeric(db$DR2)

# Filter out based on DR2 > ("threshold")
dbF <- db %>% filter(DR2 > as.numeric(args[2]))

# Group by chromosome and calculate summary statistics
summary_stats_by_chr <- dbF %>%
  group_by(chr) %>%
  summarise(
    min_DR2 = min(DR2),
    max_DR2 = max(DR2),
    mean_DR2 = mean(DR2),
    median_DR2 = median(DR2),
    min_MAF = min(MAF),
    max_MAF = max(MAF),
    mean_MAF = mean(MAF),
    median_MAF = median(MAF)
  )

overall_summary_stats <- summarise(
  dbF,
  chr = "Overall",
  min_DR2 = min(DR2),
  max_DR2 = max(DR2),
  mean_DR2 = mean(DR2),
  median_DR2 = median(DR2),
  min_MAF = min(MAF),
  max_MAF = max(MAF),
  mean_MAF = mean(MAF),
  median_MAF = median(MAF)
)

# Add overall summary statistics to the table
summary_stats_by_chr <- rbind(summary_stats_by_chr, overall_summary_stats)


# Print the summary statistics for each chromosome
write.table(summary_stats_by_chr, args[3])

options(bitmapType='cairo')
tiff(paste (c(args[4], '.tiff'), collapse=''), height =16, width = 16, units = 'cm', compression="lzw", res=600,pointsize=2,family = "Arial")
p <- ggplot(data = dbF, aes(x = MAF, y = DR2)) + geom_point()
p + facet_wrap(~chr)
dev.off()

# Generate correlation plot for the full dataset
p_full <- ggplot(data = dbF, aes(x = MAF, y = DR2)) + geom_point()
ggsave(paste0(args[4], '_full.tiff'), plot = p_full, device = "tiff", height = 8.9, width = 8.9, units = 'cm', compression = "lzw", dpi = 600)
