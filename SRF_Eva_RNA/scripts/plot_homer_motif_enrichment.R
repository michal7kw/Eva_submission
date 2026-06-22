#!/usr/bin/env Rscript

#===============================================================================
# SCRIPT: plot_homer_motif_enrichment.R
# PURPOSE: Create scatter plot showing motif enrichment vs significance
#          highlighting TEAD binding motifs
#
# DESCRIPTION:
# Creates a scatter plot similar to the reference image showing:
# - X-axis: Log P-value (significance)
# - Y-axis: % Enrichment (target % - background %)
# - TEAD motifs labeled and highlighted
# - Color gradient showing density/significance
# - Significance threshold line
#
# INPUTS:
# - HOMER knownResults.txt files from motif analysis
#
# OUTPUTS:
# - PDF and PNG plots for each direction (up/down/all)
# - Combined comparison plot
#
# USAGE:
# Rscript scripts/plot_homer_motif_enrichment.R
#===============================================================================

# Load libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrepel)
    library(dplyr)
    library(tidyr)
    library(viridis)
    library(patchwork)
})

#===============================================================================
# CONFIGURATION
#===============================================================================

BASE_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_RNA"
HOMER_DIR <- file.path(BASE_DIR, "results/07_homer_motifs")
OUTPUT_DIR <- file.path(HOMER_DIR, "plots")

# Create output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Significance threshold (for vertical line)
PVAL_THRESHOLD <- -log10(0.05)  # ~1.3 on log10 scale
LOG_PVAL_THRESHOLD <- log(0.05)  # ~-3 on natural log scale

# Motif families to highlight
HIGHLIGHT_MOTIFS <- c("TEAD", "TEA")

# Output settings
DPI <- 300

# Text size scaling factor (increase this to make all text bigger)
# Default is 1.0, use 1.5 for 50% larger text, 2.0 for double size, etc.
TEXT_SCALE <- 2.0

# Base font sizes (will be multiplied by TEXT_SCALE)
BASE_SIZE <- 12 * TEXT_SCALE        # Base font size for theme
TITLE_SIZE <- 14 * TEXT_SCALE       # Plot title size
SUBTITLE_SIZE <- 11 * TEXT_SCALE    # Plot subtitle size
LABEL_SIZE <- 3.5 * TEXT_SCALE      # TEAD label size on scatter plots
AXIS_TITLE_SIZE <- 12 * TEXT_SCALE  # Axis title size
LEGEND_SIZE <- 10 * TEXT_SCALE      # Legend text size
POINT_SIZE_TEAD <- 4 * TEXT_SCALE   # Size of TEAD points
POINT_SIZE_OTHER <- 2              # Size of other points (keep smaller)

#===============================================================================
# FUNCTIONS
#===============================================================================

read_homer_known_results <- function(filepath) {
    #' Read HOMER knownResults.txt file
    #'
    #' @param filepath Path to knownResults.txt
    #' @return Data frame with parsed results

    if (!file.exists(filepath)) {
        warning(paste("File not found:", filepath))
        return(NULL)
    }

    # Read the file
    df <- read.delim(filepath, header = TRUE, stringsAsFactors = FALSE,
                     check.names = FALSE)

    # Rename columns for easier handling
    colnames(df) <- c("motif_name", "consensus", "pvalue", "log_pvalue",
                      "qvalue", "target_count", "target_pct",
                      "bg_count", "bg_pct")

    # Clean up percentage columns (remove % sign and convert to numeric)
    df$target_pct <- as.numeric(gsub("%", "", df$target_pct))
    df$bg_pct <- as.numeric(gsub("%", "", df$bg_pct))

    # Calculate enrichment (difference in percentages)
    df$enrichment <- df$target_pct - df$bg_pct

    # Extract motif family (first part before parenthesis or slash)
    df$motif_family <- gsub("\\(.*", "", df$motif_name)
    df$motif_family <- gsub("/.*", "", df$motif_family)

    # Check if motif is TEAD-related
    df$is_tead <- grepl("TEAD|TEA\\)", df$motif_name, ignore.case = TRUE)

    # Create short label for TEAD motifs
    df$label <- ifelse(df$is_tead,
                       gsub("\\(.*", "", df$motif_name),
                       NA)

    # Add significance flag
    df$significant <- df$qvalue < 0.05

    return(df)
}


create_enrichment_plot <- function(df, title = "Motif Enrichment",
                                   show_legend = TRUE) {
    #' Create scatter plot of enrichment vs significance
    #'
    #' @param df Data frame from read_homer_known_results
    #' @param title Plot title
    #' @param show_legend Whether to show legend
    #' @return ggplot object

    if (is.null(df) || nrow(df) == 0) {
        return(NULL)
    }

    # Create significance groups for coloring
    df$sig_group <- case_when(
        df$is_tead & df$significant ~ "TEAD (Significant)",
        df$is_tead & !df$significant ~ "TEAD (Not Sig.)",
        df$significant ~ "Significant",
        TRUE ~ "Not Significant"
    )

    # Define colors
    colors <- c(
        "TEAD (Significant)" = "#E31A1C",  # Red
        "TEAD (Not Sig.)" = "#FB9A99",     # Light red
        "Significant" = "#FF7F00",          # Orange
        "Not Significant" = "#A6A6A6"       # Gray
    )

    # Create the plot
    p <- ggplot(df, aes(x = log_pvalue, y = enrichment)) +
        # Add significance threshold line
        geom_vline(xintercept = LOG_PVAL_THRESHOLD,
                   linetype = "dashed", color = "red", linewidth = 0.8) +
        # Add horizontal line at 0
        geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
        # Plot non-TEAD points first (behind)
        geom_point(data = df %>% filter(!is_tead),
                   aes(color = sig_group),
                   alpha = 0.6, size = 2) +
        # Plot TEAD points on top
        geom_point(data = df %>% filter(is_tead),
                   aes(color = sig_group),
                   alpha = 0.9, size = 3.5) +
        # Add labels for TEAD motifs
        geom_text_repel(
            data = df %>% filter(is_tead),
            aes(label = label),
            size = 3,
            fontface = "bold",
            box.padding = 0.5,
            point.padding = 0.3,
            segment.color = "gray50",
            segment.size = 0.3,
            max.overlaps = 20,
            min.segment.length = 0
        ) +
        # Color scale
        scale_color_manual(values = colors, name = "Motif Category") +
        # Labels
        labs(
            title = title,
            x = "Log P-value",
            y = "% Enrichment (Target - Background)"
        ) +
        # Add text annotations for significance regions
        # annotate("text", x = LOG_PVAL_THRESHOLD + 1, y = max(df$enrichment) * 0.95,
        #          label = "Not Sig.", color = "red", fontface = "italic", size = 3.5) +
        # annotate("text", x = LOG_PVAL_THRESHOLD - 3, y = max(df$enrichment) * 0.95,
        #          label = "Sign.", color = "red", fontface = "italic", size = 3.5) +
        # Theme
        theme_classic(base_size = 12) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            legend.position = if(show_legend) "right" else "none",
            legend.background = element_rect(fill = "white", color = "gray80"),
            panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
            axis.line = element_line(color = "black"),
            axis.ticks = element_line(color = "black")
        )

    return(p)
}


create_density_enrichment_plot <- function(df, title = "Motif Enrichment") {
    #' Create scatter plot with density coloring (like reference image)
    #'
    #' @param df Data frame from read_homer_known_results
    #' @param title Plot title
    #' @return ggplot object

    if (is.null(df) || nrow(df) == 0) {
        return(NULL)
    }

    # Calculate point density for coloring
    get_density <- function(x, y, n = 100) {
        dens <- MASS::kde2d(x, y, n = n)
        ix <- findInterval(x, dens$x)
        iy <- findInterval(y, dens$y)
        ii <- cbind(ix, iy)
        return(dens$z[ii])
    }

    # Add density (handle edge cases)
    df$density <- tryCatch({
        get_density(df$log_pvalue, df$enrichment)
    }, error = function(e) {
        rep(1, nrow(df))
    })

    # Separate TEAD and non-TEAD for layering
    df_tead <- df %>% filter(is_tead)
    df_other <- df %>% filter(!is_tead)

    # Create the plot
    p <- ggplot() +
        # Add significance threshold line
        geom_vline(xintercept = LOG_PVAL_THRESHOLD,
                   linetype = "dashed", color = "red", linewidth = 0.8) +
        # Add horizontal line at 0
        geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
        # Plot non-TEAD points with density coloring
        geom_point(data = df_other,
                   aes(x = log_pvalue, y = enrichment, color = density),
                   alpha = 0.7, size = 2) +
        # Color scale for density
        scale_color_gradientn(
            colors = c("gray70", "gray50", "orange", "red"),
            name = "Density"
        ) +
        # Add new color scale for TEAD points
        ggnewscale::new_scale_color() +
        # Plot TEAD points in distinct color
        geom_point(data = df_tead,
                   aes(x = log_pvalue, y = enrichment,
                       color = ifelse(significant, "TEAD (Sig.)", "TEAD")),
                   size = 4, alpha = 0.9) +
        scale_color_manual(
            values = c("TEAD (Sig.)" = "#1F78B4", "TEAD" = "#A6CEE3"),
            name = "TEAD Motifs"
        ) +
        # Add labels for TEAD motifs
        geom_text_repel(
            data = df_tead,
            aes(x = log_pvalue, y = enrichment, label = label),
            size = 3,
            fontface = "bold",
            color = "#1F78B4",
            box.padding = 0.5,
            point.padding = 0.3,
            segment.color = "gray50",
            segment.size = 0.3,
            max.overlaps = 20,
            min.segment.length = 0
        ) +
        # Labels
        labs(
            title = title,
            x = "Log P-value",
            y = "% Enrichment"
        ) +
        # Add text annotations
        # annotate("text", x = LOG_PVAL_THRESHOLD + 1, y = max(df$enrichment, na.rm = TRUE) * 0.95,
        #          label = "Not Sig.", color = "red", fontface = "italic", size = 3.5) +
        # annotate("text", x = LOG_PVAL_THRESHOLD - 3, y = max(df$enrichment, na.rm = TRUE) * 0.95,
        #          label = "Sign.", color = "red", fontface = "italic", size = 3.5) +
        # Theme
        theme_classic(base_size = 12) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            legend.position = "right",
            legend.background = element_rect(fill = "white", color = "gray80"),
            panel.grid.major = element_line(color = "gray90", linewidth = 0.3)
        )

    return(p)
}


create_simple_enrichment_plot <- function(df, title = "Motif Enrichment") {
    #' Create simple scatter plot matching reference style
    #'
    #' @param df Data frame from read_homer_known_results
    #' @param title Plot title
    #' @return ggplot object

    if (is.null(df) || nrow(df) == 0) {
        return(NULL)
    }

    # Separate TEAD and non-TEAD
    df_tead <- df %>% filter(is_tead)
    df_other <- df %>% filter(!is_tead)

    # Determine color by significance
    df_other$color_group <- ifelse(df_other$significant, "Significant", "Not Significant")

    # Create base plot
    p <- ggplot() +
        # Significance threshold line
        geom_vline(xintercept = LOG_PVAL_THRESHOLD,
                   linetype = "dashed", color = "red", linewidth = 0.8 * TEXT_SCALE) +
        # Zero enrichment line
        geom_hline(yintercept = 0, linetype = "dotted", color = "gray50", linewidth = 0.5) +
        # Non-TEAD points - significant ones in orange, others in gray
        geom_point(data = df_other %>% filter(!significant),
                   aes(x = log_pvalue, y = enrichment),
                   color = "gray60", alpha = 0.5, size = POINT_SIZE_OTHER) +
        geom_point(data = df_other %>% filter(significant),
                   aes(x = log_pvalue, y = enrichment),
                   color = "darkorange", alpha = 0.7, size = POINT_SIZE_OTHER) +
        # TEAD points - larger and blue
        geom_point(data = df_tead,
                   aes(x = log_pvalue, y = enrichment),
                   color = "#2166AC", size = POINT_SIZE_TEAD, alpha = 0.9) +
        # Labels for TEAD motifs
        geom_text_repel(
            data = df_tead,
            aes(x = log_pvalue, y = enrichment, label = label),
            size = LABEL_SIZE,
            fontface = "bold",
            color = "#2166AC",
            box.padding = 0.6,
            point.padding = 0.4,
            segment.color = "gray40",
            segment.size = 0.4,
            max.overlaps = 15,
            min.segment.length = 0,
            force = 2
        ) +
        # Labels
        labs(
            title = title,
            x = "Log P-value",
            y = "% of Enrichment"
        ) +
        # Theme matching reference
        theme_bw(base_size = BASE_SIZE) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = TITLE_SIZE),
            axis.title = element_text(size = AXIS_TITLE_SIZE, face = "bold"),
            axis.text = element_text(size = BASE_SIZE * 0.9),
            panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
            panel.grid.minor = element_blank(),
            legend.position = "none"
        )

    return(p)
}


#===============================================================================
# MAIN ANALYSIS
#===============================================================================

cat("============================================\n")
cat("HOMER Motif Enrichment Visualization\n")
cat("============================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Output directory:", OUTPUT_DIR, "\n\n")

# Check for ggnewscale package (optional)
has_ggnewscale <- requireNamespace("ggnewscale", quietly = TRUE)
if (!has_ggnewscale) {
    cat("Note: ggnewscale package not available, using simplified plots\n\n")
}

# Process each condition (DEGs only)
conditions <- c("downregulated_promoters", "upregulated_promoters", "all_degs_promoters")
plots_list <- list()

for (condition in conditions) {
    cat("Processing:", condition, "\n")

    # Read data
    input_file <- file.path(HOMER_DIR, condition, "knownResults.txt")
    df <- read_homer_known_results(input_file)

    if (is.null(df)) {
        cat("  - Skipping: file not found\n")
        next
    }

    cat("  - Loaded", nrow(df), "motifs\n")

    # Count TEAD motifs
    tead_count <- sum(df$is_tead)
    tead_sig <- sum(df$is_tead & df$significant)
    cat("  - TEAD motifs:", tead_count, "(", tead_sig, "significant)\n")

    # Create nice title
    nice_title <- gsub("_promoters", "", condition)
    nice_title <- gsub("_", " ", nice_title)
    nice_title <- paste0(toupper(substring(nice_title, 1, 1)),
                         substring(nice_title, 2))
    if (condition == "all_promoters") {
        nice_title <- "All Gene Promoters (Baseline)"
    } else {
        nice_title <- paste0(nice_title, " Genes - Promoter Motifs")
    }

    # Create the simple plot (matches reference style)
    p <- create_simple_enrichment_plot(df, title = nice_title)

    if (!is.null(p)) {
        plots_list[[condition]] <- p

        # Save individual plot
        pdf_file <- file.path(OUTPUT_DIR, paste0(condition, "_enrichment_plot.pdf"))
        png_file <- file.path(OUTPUT_DIR, paste0(condition, "_enrichment_plot.png"))

        ggsave(pdf_file, p, width = 10, height = 8)
        ggsave(png_file, p, width = 10, height = 8, dpi = DPI)

        cat("  - Saved:", basename(pdf_file), "\n")
    }
}

# Create combined comparison plot
if (length(plots_list) >= 2) {
    cat("\nCreating combined comparison plot...\n")

    # Combine plots in 2x2 grid
    combined <- wrap_plots(plots_list, ncol = 2) +
        plot_annotation(
            title = "HOMER Motif Enrichment Analysis - TEAD Binding Motifs",
            subtitle = "TES vs GFP: DEG Promoter Motif Analysis",
            theme = theme(
                plot.title = element_text(hjust = 0.5, face = "bold", size = TITLE_SIZE * 1.2),
                plot.subtitle = element_text(hjust = 0.5, size = SUBTITLE_SIZE)
            )
        )

    ggsave(file.path(OUTPUT_DIR, "combined_enrichment_comparison.pdf"),
           combined, width = 16, height = 14)
    ggsave(file.path(OUTPUT_DIR, "combined_enrichment_comparison.png"),
           combined, width = 16, height = 14, dpi = DPI)

    cat("  - Saved: combined_enrichment_comparison.pdf/png\n")
}

#===============================================================================
# TEAD-FOCUSED SUMMARY PLOT
#===============================================================================

cat("\nCreating TEAD-focused summary...\n")

# Combine all TEAD results for comparison
tead_summary <- data.frame()

for (condition in conditions) {
    input_file <- file.path(HOMER_DIR, condition, "knownResults.txt")
    df <- read_homer_known_results(input_file)

    if (!is.null(df)) {
        tead_df <- df %>%
            filter(is_tead) %>%
            mutate(condition = gsub("_promoters", "", condition))
        tead_summary <- rbind(tead_summary, tead_df)
    }
}

if (nrow(tead_summary) > 0) {
    # Order conditions for proper display
    tead_summary$condition <- factor(tead_summary$condition,
                                     levels = c("all_degs", "upregulated", "downregulated"))

    # Create TEAD comparison plot - bar chart
    p_tead <- ggplot(tead_summary, aes(x = condition, y = -log_pvalue, fill = condition)) +
        geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
        facet_wrap(~label, scales = "free_y", ncol = 3) +
        scale_fill_manual(values = c("downregulated" = "#2166AC", # Blue
                                     "upregulated" = "#B2182B",   # Red
                                     "all_degs" = "#762A83"),     # Purple
                          labels = c("downregulated" = "Downregulated",
                                     "upregulated" = "Upregulated",
                                     "all_degs" = "All DEGs")) +
        labs(
            title = "TEAD Motif Significance by DEG Category",
            subtitle = "Higher values indicate stronger enrichment (red line = significance threshold)",
            x = "Gene Category",
            y = "-Log P-value (Significance)",
            fill = "Category"
        ) +
        theme_minimal(base_size = BASE_SIZE) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = TITLE_SIZE),
            plot.subtitle = element_text(hjust = 0.5, size = SUBTITLE_SIZE),
            axis.title = element_text(size = AXIS_TITLE_SIZE),
            axis.text = element_text(size = BASE_SIZE * 0.9),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom",
            legend.text = element_text(size = LEGEND_SIZE),
            legend.title = element_text(size = LEGEND_SIZE),
            strip.text = element_text(face = "bold", size = BASE_SIZE)
        ) +
        geom_hline(yintercept = -LOG_PVAL_THRESHOLD, linetype = "dashed",
                   color = "red", linewidth = 0.5 * TEXT_SCALE)

    ggsave(file.path(OUTPUT_DIR, "tead_motif_comparison.pdf"),
           p_tead, width = 14, height = 8)
    ggsave(file.path(OUTPUT_DIR, "tead_motif_comparison.png"),
           p_tead, width = 14, height = 8, dpi = DPI)

    cat("  - Saved: tead_motif_comparison.pdf/png\n")

    # Create a focused TEAD dominance plot - showing enrichment difference
    tead_pivot <- tead_summary %>%
        select(label, condition, log_pvalue, enrichment) %>%
        mutate(neg_log_p = -log_pvalue)

    p_tead_dominance <- ggplot(tead_pivot, aes(x = label, y = neg_log_p, fill = condition)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.85) +
        scale_fill_manual(values = c("downregulated" = "#2166AC",
                                     "upregulated" = "#B2182B",
                                     "all_degs" = "#762A83"),
                          labels = c("downregulated" = "Downregulated DEGs",
                                     "upregulated" = "Upregulated DEGs",
                                     "all_degs" = "All DEGs")) +
        geom_hline(yintercept = -LOG_PVAL_THRESHOLD, linetype = "dashed",
                   color = "red", linewidth = 0.8 * TEXT_SCALE) +
        annotate("text", x = 0.5, y = -LOG_PVAL_THRESHOLD + 1,
                 label = "Significance threshold (p=0.05)", hjust = 0,
                 color = "red", size = LABEL_SIZE * 0.8, fontface = "italic") +
        labs(
            title = "TEAD Motif Enrichment: Dominance in Downregulated Genes",
            subtitle = "TEAD motifs are specifically enriched in TES-repressed genes",
            x = "TEAD Motif",
            y = "-Log P-value (Significance)",
            fill = "DEG Category"
        ) +
        theme_minimal(base_size = BASE_SIZE) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = TITLE_SIZE),
            plot.subtitle = element_text(hjust = 0.5, size = SUBTITLE_SIZE),
            axis.title = element_text(size = AXIS_TITLE_SIZE),
            axis.text = element_text(size = BASE_SIZE * 0.9),
            axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
            legend.position = "right",
            legend.text = element_text(size = LEGEND_SIZE),
            legend.title = element_text(size = LEGEND_SIZE),
            panel.grid.minor = element_blank()
        ) +
        coord_flip()

    ggsave(file.path(OUTPUT_DIR, "tead_dominance_plot.pdf"),
           p_tead_dominance, width = 12, height = 8)
    ggsave(file.path(OUTPUT_DIR, "tead_dominance_plot.png"),
           p_tead_dominance, width = 12, height = 8, dpi = DPI)

    cat("  - Saved: tead_dominance_plot.pdf/png\n")

    # Save TEAD summary table
    tead_table <- tead_summary %>%
        select(condition, motif_name, consensus, log_pvalue, qvalue,
               target_pct, bg_pct, enrichment) %>%
        arrange(condition, log_pvalue)

    write.csv(tead_table, file.path(OUTPUT_DIR, "tead_motif_summary.csv"),
              row.names = FALSE)
    cat("  - Saved: tead_motif_summary.csv\n")
}

#===============================================================================
# COMPLETION
#===============================================================================

cat("\n============================================\n")
cat("Analysis Complete!\n")
cat("============================================\n")
cat("\nOutput files saved to:", OUTPUT_DIR, "\n")
cat("\nKey findings:\n")

# Print TEAD significance summary
for (condition in conditions) {
    input_file <- file.path(HOMER_DIR, condition, "knownResults.txt")
    df <- read_homer_known_results(input_file)

    if (!is.null(df)) {
        tead_sig <- df %>% filter(is_tead & significant)
        nice_name <- gsub("_promoters", "", condition)

        if (nrow(tead_sig) > 0) {
            cat("\n", nice_name, ":\n", sep = "")
            cat("  TEAD motifs significantly enriched:\n")
            for (i in 1:min(5, nrow(tead_sig))) {
                cat("    -", tead_sig$label[i],
                    "(enrichment:", round(tead_sig$enrichment[i], 2),
                    "%, p=", format(tead_sig$pvalue[i], scientific = TRUE), ")\n")
            }
        } else {
            cat("\n", nice_name, ": No significant TEAD enrichment\n", sep = "")
        }
    }
}

cat("\n")
