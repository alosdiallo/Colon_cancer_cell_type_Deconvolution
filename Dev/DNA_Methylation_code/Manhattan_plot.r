```{r}

library(qqman)
# Create a vector of indices to match rownames of fit2$p.value with the rows of annotation
matching_indices <- match(rownames(fit2$p.value), annotation$probeID)

# Extract the corresponding chromosome and position data
chromosome_data <- annotation$CpG_chrm[matching_indices]
position_data <- annotation$CpG_beg[matching_indices]

# Extract the corresponding p-values
p_values <- fit2$p.value[, "Post_Drug - Baseline"]

# Combine into a dataframe
cpg_data <- data.frame(
  chromosome = chromosome_data,
  position = position_data,
  p.value = p_values
)

# Check the first few rows of cpg_data to ensure it looks correct
head(cpg_data)

# Remove 'chr' prefix from the chromosome column
cpg_data$chromosome <- gsub("chr", "", cpg_data$chromosome)

# Convert 'X', 'Y', and 'MT' to numeric values
cpg_data$chromosome[cpg_data$chromosome == "X"] <- "23"
cpg_data$chromosome[cpg_data$chromosome == "Y"] <- "24"
cpg_data$chromosome[cpg_data$chromosome == "MT"] <- "25"

# Convert chromosome to numeric
cpg_data$chromosome <- as.numeric(cpg_data$chromosome)

# Check for NA values again after the conversion
sum(is.na(cpg_data$chromosome))

# Calculate -log10(p-value) for the Manhattan plot
cpg_data$minuslog10p <- -log10(cpg_data$p.value)

threshold <- -log10(5e-2)  # The genome-wide significance line
significant_points <- cpg_data[cpg_data$minuslog10p > threshold, ]

# Remove any rows with NA or infinite p-values
cpg_data <- cpg_data[!is.na(cpg_data$minuslog10p) & !is.infinite(cpg_data$minuslog10p), ]
cpg_data$SNP <- rownames(cpg_data)
# Now create the Manhattan plot
# Assuming no NA values, run the manhattan function again
manhattan(cpg_data,
          chr = "chromosome",
          bp = "position",
          p = "minuslog10p",
          col = c("blue4", "orange3"),
          main = "Manhattan Plot for Differential Methylation",
          ylim = c(0, max(cpg_data$minuslog10p, na.rm = TRUE)),
          suggestiveline = -log10(1e-4),
          genomewideline = -log10(5e-2),
          annotatePval = 0.01
)
#with(significant_points, text(position, minuslog10p, labels=rownames(significant_points), pos=3, cex=0.6))
