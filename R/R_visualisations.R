# This script generates plots for Crn4 analysis including tree alignment plots,
# fusion plots, and distribution tables

# ==============================================================================
# 1. SETUP AND CONFIGURATION
# ==============================================================================

project <- "211024" # Create a folder with this name in the projects directory

# Install and load required libraries.
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("ggtree", quietly = TRUE)) install.packages("ggtree")
if (!requireNamespace("gggenes", quietly = TRUE)) install.packages("gggenes")
if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
if (!requireNamespace("phangorn", quietly = TRUE)) install.packages("phangorn")
if (!requireNamespace("ape", quietly = TRUE)) install.packages("ape")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(ggplot2)
library(ggtree)
library(gggenes)
library(gridExtra)
library(patchwork)
library(phangorn)
library(ape)
library(dplyr)

# ==============================================================================
# 2. DATA LOADING
# ==============================================================================

# Load tree data
crn4_tree_path <- paste("projects", project, "data", "crn4_tree.txt", sep = "/")
crn4_tree <- read.tree(crn4_tree_path)

# Edit the tip labels to remove everything after __
crn4_tree$tip.label <- gsub("__.*", "", crn4_tree$tip.label)
# Then remove everything after second occurrence of underscore
crn4_tree$tip.label <- gsub("_[^_]*$", "", crn4_tree$tip.label)

# Remove specific nodes from tree
crn4_tree <- drop.tip(crn4_tree, "GCF_002285795.1")
crn4_tree <- drop.tip(crn4_tree, "GCF_002863905.1")

# Root the tree
crn4_tree_rooted <- root(crn4_tree, outgroup = "GCF_009930875.1", resolve.root = TRUE)
rooted_tree_final <- crn4_tree_rooted  # This variable is used later in the plotting

# Load alignment data
alignment_path <- paste("projects", project, "data", "crn4.afa", sep = "/")
alignment <- read.FASTA(alignment_path, type = "AA")

# Convert to phydat object
alignment_phydat <- phyDat(as.character(alignment), type = "AA")

# Fit the tree to the alignment
fit_k1 <- pml(crn4_tree, alignment_phydat)  # Fit the tree to the alignment
fit_k4 <- pml(crn4_tree, alignment_phydat, k = 4)  # Fit with 4 rate categories
fit_k8 <- pml(crn4_tree, alignment_phydat, k = 8)  # Fit with 8 rate categories

# Load and process GFF data
gff_dir <- paste("projects", project, "data/gff_10k/all", sep = "/")
crn4_genes <- gff_to_gggenes(gff_dir)
crn4_genes <- subset(crn4_genes, type == "CDS")
colnames(crn4_genes)[2] <- "molecule"

# Set columns start and end to numeric
crn4_genes$start <- as.numeric(crn4_genes$start)
crn4_genes$end <- as.numeric(crn4_genes$end)
crn4_genes <- na.omit(crn4_genes)

# Create column orientation based on direction
crn4_genes$orientation <- ifelse(crn4_genes$direction == 1, 1, 0)

# Only retain the following columns
crn4_genes <- crn4_genes[, c("molecule", "gene", "start", "end", "strand", "orientation")]

# Remove everything after last underscore in column molecule
crn4_genes$molecule <- gsub("_[^_]*$", "", crn4_genes$molecule)

# ==============================================================================
# 3. GENE DATA PROCESSING
# ==============================================================================

crn4_gene_simplified <- crn4_genes

# Simplify gene names
crn4_gene_simplified$gene <- ifelse(grepl("crn4a", crn4_gene_simplified$gene), "crn4", crn4_gene_simplified$gene)
crn4_gene_simplified$gene <- ifelse(grepl("crn4b", crn4_gene_simplified$gene), "crn4", crn4_gene_simplified$gene)
crn4_gene_simplified$gene <- ifelse(grepl("Cas10", crn4_gene_simplified$gene), "Cas10", crn4_gene_simplified$gene)
crn4_gene_simplified$gene <- ifelse(grepl("Cas1", crn4_gene_simplified$gene), "Cas1", crn4_gene_simplified$gene)
crn4_gene_simplified$gene <- ifelse(grepl("Cas2", crn4_gene_simplified$gene), "Cas2", crn4_gene_simplified$gene)
crn4_gene_simplified$gene <- ifelse(grepl("Cas6", crn4_gene_simplified$gene), "Cas6", crn4_gene_simplified$gene)
crn4_gene_simplified$gene <- ifelse(grepl("can3", crn4_gene_simplified$gene), "Csm6-2", crn4_gene_simplified$gene)
crn4_gene_simplified$gene <- ifelse(grepl("cora", crn4_gene_simplified$gene), "CorA", crn4_gene_simplified$gene)
crn4_gene_simplified$gene <- ifelse(grepl("Cmr", crn4_gene_simplified$gene), "Cmr", crn4_gene_simplified$gene)
crn4_gene_simplified$gene <- ifelse(grepl("crn4", crn4_gene_simplified$gene), "Crn4", crn4_gene_simplified$gene)
crn4_gene_simplified$gene <- ifelse(grepl("Csm,Other", crn4_gene_simplified$gene), "Csm", crn4_gene_simplified$gene)
crn4_gene_simplified$gene <- ifelse(grepl("can1-2", crn4_gene_simplified$gene), "Can2", crn4_gene_simplified$gene)
crn4_gene_simplified$gene <- ifelse(grepl("Other", crn4_gene_simplified$gene), "Other", crn4_gene_simplified$gene)

# Create category column
crn4_gene_simplified$category <- ifelse(grepl("Cas", crn4_gene_simplified$gene), "Cas", crn4_gene_simplified$gene)
crn4_gene_simplified$category <- ifelse(grepl("Csm", crn4_gene_simplified$gene), "Cas", crn4_gene_simplified$category)
crn4_gene_simplified$category <- ifelse(grepl("Cmr", crn4_gene_simplified$gene), "Cas", crn4_gene_simplified$category)
crn4_gene_simplified$category <- ifelse(grepl("Crn4", crn4_gene_simplified$gene), "Crn4", crn4_gene_simplified$category)
crn4_gene_simplified$category <- ifelse(grepl("Csm6-2", crn4_gene_simplified$gene), "Csm6-2", crn4_gene_simplified$category)
crn4_gene_simplified$category <- ifelse(grepl("CorA", crn4_gene_simplified$gene), "CorA", crn4_gene_simplified$category)
crn4_gene_simplified$category <- ifelse(grepl("Can2", crn4_gene_simplified$gene), "Can2", crn4_gene_simplified$category)

# Create labels
crn4_gene_simplified$label <- crn4_gene_simplified$category
crn4_gene_simplified$label <- as.character(crn4_gene_simplified$label)
crn4_gene_simplified$label <- ifelse(grepl("Other", crn4_gene_simplified$label), "", crn4_gene_simplified$label)

# Handle specific fusion proteins
crn4_gene_simplified$category <- ifelse(crn4_gene_simplified$molecule == "GCF_003111765.1" & 
                                        crn4_gene_simplified$gene == "Csm6-2" & 
                                        crn4_gene_simplified$start == 517169, 
                                        "CorA/Csm6-2 fusion", crn4_gene_simplified$category)
crn4_gene_simplified$category <- ifelse(crn4_gene_simplified$molecule == "GCF_004102025.1" & 
                                        crn4_gene_simplified$gene == "Csm6-2" & 
                                        crn4_gene_simplified$start == 520910, 
                                        "CorA/Csm6-2 fusion", crn4_gene_simplified$category)

# If orientation is 0, swap start and end
crn4_gene_simplified$temp_start <- ifelse(crn4_gene_simplified$orientation == 0, 
                                           crn4_gene_simplified$end, 
                                           crn4_gene_simplified$start)
crn4_gene_simplified$temp_end <- ifelse(crn4_gene_simplified$orientation == 0, 
                                         crn4_gene_simplified$start, 
                                         crn4_gene_simplified$end)
crn4_gene_simplified$start <- crn4_gene_simplified$temp_start
crn4_gene_simplified$end <- crn4_gene_simplified$temp_end
crn4_gene_simplified <- subset(crn4_gene_simplified, select = -c(temp_start, temp_end))

# ==============================================================================
# 4. LOAD METADATA
# ==============================================================================

# Load metadata
info <- read.csv(paste("projects", project, "data", "mastertable_v2.tsv", sep="/"), sep="\t")
info_crn4 <- subset(info, crn4a == "True" | crn4b == "True")
rownames(info_crn4) <- info_crn4$sample

# Convert cyclase to boolean
info_crn4$cyclase <- as.logical(info_crn4$cyclase)
info_crn4 <- info_crn4[, c('sample', setdiff(names(info_crn4), 'sample'))]

# ==============================================================================
# 5. PREPARE FUSION DATA
# ==============================================================================

# Create fusion dataframe
crn4_fusions <- data.frame(
  molecule = character(),
  gene = character(),
  start = numeric(),
  end = numeric(),
  strand = character(),
  orientation = numeric(),
  stringsAsFactors = FALSE
)

# Add fusion entries
crn4_fusions <- rbind(crn4_fusions, data.frame(
  molecule = "WP_021796948.1",
  gene = "Csm6-2/Crn4a",
  start = 1,
  end = 987,
  strand = "forward",
  orientation = 1
))

crn4_fusions <- rbind(crn4_fusions, data.frame(
  molecule = "WP_197719274.1",
  gene = "Csm6-2/Crn4a",
  start = 1,
  end = 957,
  strand = "forward",
  orientation = 1
))

crn4_fusions <- rbind(crn4_fusions, data.frame(
  molecule = "WP_126376302.1",
  gene = "TIR-SAVED/Crn4b",
  start = 1,
  end = 604,
  strand = "forward",
  orientation = 1
))

# Create subgene coordinates dataframe
crn4_fusion_coordinates <- data.frame(
  molecule = character(),
  gene = character(),
  start = numeric(),
  end = numeric(),
  strand = character(),
  subgene = character(),
  from = numeric(),
  to = numeric(),
  orientation = numeric(),
  stringsAsFactors = FALSE
)

# Add subgene entries for WP_021796948.1
crn4_fusion_coordinates <- rbind(crn4_fusion_coordinates, data.frame(
  molecule = "WP_021796948.1",
  gene = "Csm6-2/Crn4a",
  start = 1,
  end = 987,
  strand = "forward",
  subgene = "Crn4a",
  from = 838,
  to = 987,
  orientation = 1
))

crn4_fusion_coordinates <- rbind(crn4_fusion_coordinates, data.frame(
  molecule = "WP_021796948.1",
  gene = "Csm6-2/Crn4a",
  start = 1,
  end = 987,
  strand = "forward",
  subgene = "Csm6-2",
  from = 1,
  to = 824,
  orientation = 1
))

# Add subgene entries for WP_197719274.1
crn4_fusion_coordinates <- rbind(crn4_fusion_coordinates, data.frame(
  molecule = "WP_197719274.1",
  gene = "Csm6-2/Crn4a",
  start = 1,
  end = 957,
  strand = "forward",
  subgene = "Crn4a",
  from = 822,
  to = 957,
  orientation = 1
))

crn4_fusion_coordinates <- rbind(crn4_fusion_coordinates, data.frame(
  molecule = "WP_197719274.1",
  gene = "Csm6-2/Crn4a",
  start = 1,
  end = 957,
  strand = "forward",
  subgene = "Csm6-2",
  from = 1,
  to = 808,
  orientation = 1
))

# Add subgene entries for WP_126376302.1
crn4_fusion_coordinates <- rbind(crn4_fusion_coordinates, data.frame(
  molecule = "WP_126376302.1",
  gene = "TIR-SAVED/Crn4b",
  start = 1,
  end = 604,
  strand = "forward",
  subgene = "Crn4b",
  from = 474,
  to = 604,
  orientation = 1
))

crn4_fusion_coordinates <- rbind(crn4_fusion_coordinates, data.frame(
  molecule = "WP_126376302.1",
  gene = "TIR-SAVED/Crn4b",
  start = 1,
  end = 604,
  strand = "forward",
  subgene = "TIR-SAVED",
  from = 1,
  to = 460,
  orientation = 1
))

# Set factor levels for plotting order
crn4_fusions$molecule <- factor(crn4_fusions$molecule, 
                                levels = c("WP_021796948.1", "WP_197719274.1", "WP_126376302.1"))
crn4_fusion_coordinates$molecule <- factor(crn4_fusion_coordinates$molecule, 
                                           levels = c("WP_021796948.1", "WP_197719274.1", "WP_126376302.1"))

# Add midpoint for labeling
crn4_fusions$midpoint <- (crn4_fusions$start + crn4_fusions$end) / 2

# ==============================================================================
# 6. CREATE PLOTS
# ==============================================================================

# Set factor levels for gene categories
crn4_gene_simplified$category <- as.factor(crn4_gene_simplified$category)
crn4_gene_simplified$category <- factor(crn4_gene_simplified$category, 
                                        levels = c("Crn4", "Csm6-2", "CorA", "CorA/Csm6-2 fusion", 
                                                   "Can2", "Cas", "Other"))

# Create tree plot
crn4_tree_plot <- ggtree(rooted_tree_final)
crn4_tree_plot <- crn4_tree_plot %<+% info_crn4
crn4_tree_plot <- crn4_tree_plot + 
  geom_tiplab(align=TRUE, linetype='dashed', linesize=0.5, alpha=0.3, offset = 0.6, size = 3) +
  geom_point2(aes(subset = cyclase == FALSE), color = "red", size = 3) +
  geom_label2(aes(subset = crn4a == "True"), label = "crn4a", fill = "#E09F3E", 
              color = "black", size = 3, hjust = -0.2, vjust = 0.5) +
  geom_label2(aes(subset = crn4b == "True"), label = "crn4b", fill = "#9c6c25", 
              color = "#fff", size = 3, hjust = -0.2, vjust = 0.5) +
  geom_point2(aes(subset = locus == "GCF_004798665.1_0"), fill = "#f5bd2f", 
              color = "black", size = 3, shape = 22) +
  xlim_tree(5.5)

# Create tree alignment plot
crn4_tree_align_plot <- crn4_tree_plot +
  geom_facet(data = crn4_gene_simplified, 
             mapping = aes(xmin = start, xmax = end, fill = category, color = category),
             geom = geom_motif, panel = 'Alignment',
             label = 'category', align = 'center', on = 'Crn4') +
  scale_x_continuous(expand=c(0,0)) +
  theme(strip.text=element_blank(),
        panel.spacing=unit(0, 'cm')) +
  scale_fill_manual(values = c("Cas" = "#dbdbdb", "Crn4" = "#E09F3E", "Csm6-2" = "#c1ed98", 
                               "CorA" = "#fa9394", "CorA/Csm6-2 fusion" = "#9dba82", 
                               "Can2" = "#a6c0e3", "Other" = "#ffffff")) +
  scale_color_manual(values = c("Cas" = "#969696", "Crn4" = "#3b3b3b", "Csm6-2" = "#3b3b3b", 
                                "CorA" = "#3b3b3b", "CorA/Csm6-2 fusion" = "#000", 
                                "Can2" = "#3b3b3b", "Other" = "#969696")) +
  labs(fill = "Gene category", color = "Gene category")

# Adjust facet widths
crn4_tree_align_plot <- facet_widths(crn4_tree_align_plot, widths = c(8, 9))

# Create fusion plot
fusion_plot <- ggplot(crn4_fusions, aes(xmin = start, xmax = end, y = molecule)) +
  geom_gene_arrow(fill = "white") +
  geom_subgene_arrow(data = crn4_fusion_coordinates,
                     aes(xmin = start, xmax = end, fill = subgene,
                         xsubmin = from, xsubmax = to), color="black", alpha=.7) +
  geom_text(aes(x = midpoint, y = molecule, label = molecule),
            vjust = 0.5, size = 2.5, color = "#404040") +
  theme_genes() +
  xlim(0, 1000) +
  scale_y_discrete(limits = rev(levels(crn4_fusions$molecule))) +
  theme_genes() %+replace% theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank()
  ) +
  scale_fill_manual(values = c("TIR-SAVED" = "#c9aae6", "Crn4a" = "#E09F3E", 
                               "Crn4b" = "#9c6c25", "Csm6-2" = "#c1ed98")) +
  labs(fill = "Domain")

# Create distribution table
crn4_table <- data.frame(
  Prokaryotes = c(12, 10), 
  Phages = c(0, 11), 
  Plasmids = c(0, 10)
)
rownames(crn4_table) <- c("Crn4a", "Crn4b")

crn4_table_grob <- tableGrob(crn4_table, theme = ttheme_minimal(
  core = list(fg_params = list(fontsize = 11)),
  colhead = list(fg_params = list(fontface = "bold", fontsize = 10)),
  padding = unit(c(4, 4), "mm"),
  rowhead = list(fg_params = list(fontface = "bold", fontsize = 10))
))

# Combine all plots into Figure 1
fig1 <- crn4_tree_align_plot / (fusion_plot + crn4_table_grob) +
  plot_layout(heights = c(7, 1), widths = list(c(1), c(1, 1))) +
  plot_annotation(tag_levels = "A",
                  tag_suffix = "",
                  tag_sep = "",
                  theme = theme(plot.tag = element_text(size = 16, face = "bold")))

# ==============================================================================
# 7. SAVE PLOTS
# ==============================================================================

# Save combined figure
fig1_path <- paste("projects", project, "plots_reissued/fig1_10k.png", sep = "/")
ggsave(fig1_path, plot = fig1, bg = "white", height = 9, width = 10)

# Save individual panels to 2025 folder
crn4_tree_align_plot_path <- paste("projects", project, "plots_reissued/2025/crn4_tree_align_plot.pdf", sep = "/")
ggsave(crn4_tree_align_plot_path, plot = crn4_tree_align_plot, bg = "white", height = 8, width = 11)

fusion_plot_path <- paste("projects", project, "plots_reissued/2025/fusion_plot.pdf", sep = "/")
ggsave(fusion_plot_path, plot = fusion_plot, bg = "white", height = 3, width = 5)

crn4_table_grob_path <- paste("projects", project, "plots_reissued/2025/crn4_table_grob.pdf", sep = "/")
ggsave(crn4_table_grob_path, plot = crn4_table_grob, bg = "white", height = 3, width = 5)

# Print completion message
cat("Plot generation complete!\n")
cat("Combined figure saved to:", fig1_path, "\n")
cat("Individual plots saved to: projects/", project, "/plots_reissued/2025/\n", sep="")
