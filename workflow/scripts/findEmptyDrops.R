#Rscript.exe findEmptyDrops.R --raw_feature_bc_matrix sample/raw_feature_bc_matrix/ --sample_id Cab39L_KO --outdir .

library(optparse)
library(DropletUtils)
library(tidyverse)

option_list <- list(
  make_option(c("--raw_feature_bc_matrix"), type = "character", help = ""),
  make_option(c("--sample_id"), type = "character", help = ""),
  make_option(c("--outdir"), type = "character", help = " "),
  make_option(c("--niters"), type = "numeric", default = 10000,)
)

parseobj <- OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]")
opt <- parse_args(parseobj)
# options(bitmapType='cairo', scipen=0)
raw_feature_bc_matrix <- opt$raw_feature_bc_matrix
sample_id <- opt$sample_id
outdir <- opt$outdir

# options <- commandArgs(trailingOnly = TRUE)
# path <- "../data_raw/Cab39L_KO/raw_feature_bc_matrix/"
# path <- options[1]

sce <- read10xCounts(raw_feature_bc_matrix, col.names = TRUE)

br.out <- barcodeRanks(counts(sce))
br.out.df <- as.data.frame(br.out)
br.out.df$barcode <- colData(sce)$Barcode

cells_knee <- filter(br.out.df, total > br.out@metadata$knee) %>% nrow
cells_inflection <- filter(br.out.df, total > br.out@metadata$inflection) %>% nrow

set.seed(100)
e.out <- emptyDrops(counts(sce), niters = opt$niters)
e.out.df <- e.out %>%
  as.data.frame %>%
  mutate(
    barcode = colData(sce)$Barcode,
    is.cell = FDR <= 0.01 & !is.na(FDR))

possible_cells <- filter(e.out.df, Limited == TRUE, is.cell == FALSE)
msg <- ""
if (nrow(possible_cells) > 0) {
  msg <- sprintf("More cell in sample %s possible to detect! 
                 Consider increasing the niters parameter", sample_id)
  message(msg)
}

table <- left_join(e.out.df, br.out.df, by = "barcode") %>%
  mutate(barcode = str_sub(barcode, 1, 16)) %>%
  select(barcode, total, rank, LogProb, is.cell) %>%
  filter(total > 100)

results <- lst(
  table,
  knee_point = br.out@metadata$knee,
  inflection_point = br.out@metadata$inflection,
  cells_above_knee = cells_knee,
  cells_above_inflecion = cells_inflection,
  msg
)

plot1 <- function(data) {
  data$table %>%
    filter(!is.na(LogProb)) %>%
    ggplot(aes(total, -LogProb, color = is.cell)) +
    geom_point() +
    geom_vline(xintercept = data$knee_point, linetype = 2, colour = "dodgerblue") +
    geom_vline(xintercept = data$inflection_point, linetype = 2, colour = "forestgreen") +
    labs(title = data$sample_name)
}
emptyDrops_plot <- plot1(results)

results$emptyDrops_plot <- emptyDrops_plot

files <- c(".Rds", ".tsv", ".png") %>%
  str_c(sample_id, .) %>%
  str_c(outdir, ., sep = "/") %>%
  str_replace_all("//", "/") %>%
  set_names(c("rds", "tsv", "png")) %>%
  as.list()

c(rds = str_c(sample_id, "Rds"))

saveRDS(results, file = files$rds)
write_tsv(table, file = files$tsv)
ggsave(filename = files$png, plot = emptyDrops_plot)