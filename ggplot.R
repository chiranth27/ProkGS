library(ggplot2)

# File path
file_path <- "/home/eukprogs/Results/stream.txt"
lines <- readLines(file_path)

# Extract gene descriptions (lines 1 to 5)
gene_labels <- sub("^\\w: (.+?) \\(+\\)$", "\\1(+)", lines[1:5])

gene_letters <- substr(lines[1:5], 1, 1)

# Extract gene positions
start_line <- grep("^start_values=", lines)
positions <- as.numeric(unlist(strsplit(sub("start_values=", "", lines[start_line]), ",")))

# Extract metadata
genome_line <- grep("^Genome_Accession:", lines)
organism_line <- grep("^Organism:", lines)
pdb_line <- grep("^PDB:", lines)

genome_accession <- sub("Genome_Accession: ", "", lines[genome_line])
organism <- sub("Organism: ", "", lines[organism_line])
pdb <- sub("PDB: ", "", lines[pdb_line])

# Gene dataframe
df <- data.frame(
  Gene = gene_letters,
  Position = positions,
  X = seq(-1, 3, length.out = length(gene_letters)),  # Adjust X-coordinates for visualization
  Color = ifelse(gene_letters == "C", "red", "cyan")
)

# Annotation labels (to left of gene map)
annotations <- data.frame(
  x = rep(-1.9, length(gene_labels)),
  y = seq(2, 1, length.out = length(gene_labels)),
  label = gene_labels  # Use only gene_labels
)

# Metadata info
metadata <- data.frame(
  x = 1.3,
  y = c(-1.5, -1.8, -2.1),
  label = c(
    paste("Genome Accession:", genome_accession),
    paste("Organism:", organism),
    paste("PDB:", pdb)
  )
)

# Color legend
# Color legend
legend_box <- data.frame(
  x = c(-1.5, -1.5),
  y = c(-2.8, -3.1),
  label = c("General Proteins", "Annotated Protein"),
  color = c("cyan", "red")
)


# Plotting
p <- ggplot(df, aes(x = X, y = 0)) +
  geom_segment(data = data.frame(x = min(df$X) - 0.5, xend = max(df$X) + 0.5), 
               aes(x = x, xend = xend, y = 0, yend = 0), 
               color = "white", linewidth = 1) +
  geom_point(aes(color = Color), size = 6, shape = 21, fill = df$Color, stroke = 0.1) +
  geom_text(aes(label = Gene, y = 0.18), color = "yellow", fontface = "bold", size = 6, vjust = 0) +
  geom_text(aes(label = Position, y = -0.15), color = "white", size = 3.5, angle = 0, vjust = 1) +
  annotate("text", x = min(df$X) - 0.3, y = 0.35, label = "Upstream", 
           color = "white", fontface = "bold", size = 4, hjust = 1) +
  annotate("text", x = max(df$X) + 0.05, y = 0.35, label = "Downstream", 
           color = "white", fontface = "bold", size = 3.5, hjust = 0) +
  geom_text(data = annotations, aes(x = x, y = y, label = label), 
            color = ifelse(df$Color == "red", "red", "yellow"), fontface = "bold", size = 4.5, hjust = 0) +
  geom_text(data = metadata, aes(x = x, y = y, label = label), 
            color = c("yellow", "yellow", "red"), fontface = "bold", size = 5, hjust = 0.5) +
  geom_point(data = legend_box, aes(x = x - 0.3, y = y, color = color), 
             size = 6, shape = 21, stroke = 1.2, fill = legend_box$color) +
  geom_text(data = legend_box, aes(x = x, y = y, label = label),
            color = "yellow", fontface = "bold", size = 5, hjust = 0) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "black", color = NA),
    plot.background = element_rect(fill = "black", color = NA),
    legend.position = "none"
  ) +
  ggtitle("Genomic Annotation of PDB Entry") +
  theme(plot.title = element_text(color = "yellow", size = 22, face = "bold", hjust = 0.7))

print(p)
ggsave("/home/eukprogs/Results/genomic_annotation_plot.png", width = 720/72, height = 412/72, dpi = 72)

