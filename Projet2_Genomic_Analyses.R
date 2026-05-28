############################################################
# BING-F432 - Spatial and molecular epidemiology (ULB, 2025-2026)
# Project 2 : Genomic analyses
# Team 10  : Biefnot Alexis & Thiels Hannah
# Dataset  : DENV simulated outbreak dataset 2-3
# Purpose  : preliminary epidemiological, spatial and molecular analyses
############################################################
#
# References used in this script
# ------------------------------------------------------------
# Methods and tools:
#   - EpiEstim (Rt estimation):
#       Cori A., Ferguson N.M., Fraser C., Cauchemez S. 2013.
#       A new framework and software to estimate time-varying
#       reproduction numbers during epidemics.
#       American Journal of Epidemiology 178:1505-1512.
#       https://doi.org/10.1093/aje/kwt133
#   - Phi test for recombination:
#       Bruen T.C., Philippe H., Bryant D. 2006.
#       A simple and robust statistical test for detecting the
#       presence of recombination.
#       Genetics 172:2665-2681.
#       https://doi.org/10.1534/genetics.105.048975
#   - IQ-TREE 2/3 (maximum-likelihood phylogeny):
#       Minh B.Q., Schmidt H.A., Chernomor O., Schrempf D.,
#       Woodhams M.D., von Haeseler A., Lanfear R. 2020.
#       IQ-TREE 2: new models and efficient methods for
#       phylogenetic inference in the genomic era.
#       Molecular Biology and Evolution 37:1530-1534.
#       https://doi.org/10.1093/molbev/msaa015
#   - TempEst (temporal signal):
#       Rambaut A., Lam T.T., Carvalho L.M., Pybus O.G. 2016.
#       Exploring the temporal structure of heterochronous
#       sequences using TempEst.
#       Virus Evolution 2:vew007.
#       https://doi.org/10.1093/ve/vew007
#   - BEAST X (Bayesian phylogenetic inference):
#       Suchard M.A., Lemey P., Baele G., Ayres D.L., Drummond A.J.,
#       Rambaut A. 2018.
#       Bayesian phylogenetic and phylodynamic data integration
#       using BEAST 1.10.
#       Virus Evolution 4:vey016.
#       https://doi.org/10.1093/ve/vey016
#   - Discrete phylogeographic inference:
#       Lemey P., Rambaut A., Drummond A.J., Suchard M.A. 2009.
#       Bayesian phylogeography finds its roots.
#       PLoS Computational Biology 5:e1000520.
#       https://doi.org/10.1371/journal.pcbi.1000520
#   - Skygrid coalescent model:
#       Gill M.S., Lemey P., Faria N.R., Rambaut A., Shapiro B.,
#       Suchard M.A. 2013.
#       Improving Bayesian population dynamics inference: a
#       coalescent-based model for multiple loci.
#       Molecular Biology and Evolution 30:713-724.
#       https://doi.org/10.1093/molbev/mss265
#   - Sampling bias in phylogeography (key limitation to discuss):
#       Lemey P., Rambaut A., Bedford T., Faria N., Bielejec F.,
#       Baele G., Russell C.A., Smith D.J., Pybus O.G.,
#       Brockmann D., Suchard M.A. 2014.
#       Unifying viral genetics and human transportation data
#       to predict the global transmission dynamics of human
#       influenza H3N2.
#       PLoS Pathogens 10:e1003932.
#       https://doi.org/10.1371/journal.ppat.1003932
#
# Dengue-specific epidemiological parameters:
#   - Serial interval distribution:
#       ten Bosch Q.A., Clapham H.E., Lambrechts L., Duong V.,
#       Buchy P., Althouse B.M., Lloyd A.L., Waller L.A.,
#       Morrison A.C., Kitron U., Vazquez-Prokopec G.M.,
#       Scott T.W., Perkins T.A. 2018.
#       Contributions from the silent majority dominate
#       dengue virus transmission.
#       PLoS Pathogens 14:e1006965.
#       https://doi.org/10.1371/journal.ppat.1006965
#
# R packages:
#   - ape:   Paradis E., Schliep K. 2019. Bioinformatics 35:526-528.
#   - sf:    Pebesma E. 2018. R Journal 10:439-446.
#
# Inspiration for figure design:
#   - Candido D.S. et al. 2020. Evolution and epidemic spread of
#       SARS-CoV-2 in Brazil. Science 369:1255-1260
#       (cumulative case curves, sampling lag dumbbell plot,
#        sequencing-effort scatter, choropleth maps).
#   - Faria N.R. et al. 2018. Genomic and epidemiological
#       monitoring of yellow fever virus transmission potential.
#       Science 361:894-899
#       (spatial distribution panels, time series of cases by host).
#
############################################################

# ==============================================================================
# 0. PROJECT SETUP
# ==============================================================================


# Run this script from the project root directory.
# The project root must contain data/ and outputs/.
# Prrint the current working directory
cat("Working directory:", getwd(), "\n")

# Stop immediately if the expected project structure is not present.
if (!dir.exists("data/sequences")) {
  stop("Folder not found: data/sequences. Run the script from the project root folder.")
}
if (!dir.exists("data/shapefiles")) {
  stop("Folder not found: data/shapefiles. Run the script from the project root folder.")
}

# Create the output directory if it does not already exist
dir.create("outputs", showWarnings = FALSE, recursive = TRUE)

# Required R packages for the analysis
required_packages <- c(
  "ape", "phangorn", "raster", "sf", "sp", "RColorBrewer",
  "colorspace", "lubridate", "EpiEstim", "seqinr"
)

# Install missing packages, then load all required packages
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

# ==============================================================================
# 0.1 HELPER FUNCTIONS
# ==============================================================================

# Check that a file exists before trying to read it.
check_file_exists <- function(path) {
  if (!file.exists(path)) {
    stop("Required file not found: ", path)
  }
}

# Check that expected columns are present in a data frame.
check_required_columns <- function(data, required_cols, data_name) {
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required column(s) in ", data_name, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }
}

# Convert character dates safely to Date format.
convert_dates_safely <- function(date_vector, column_name) {
  converted_dates <- lubridate::ymd(date_vector)
  
  if (any(is.na(converted_dates))) {
    bad_values <- unique(date_vector[is.na(converted_dates)])
    stop(
      "Some values in ", column_name, " could not be converted to Date format: ",
      paste(bad_values, collapse = ", ")
    )
  }
  
  return(converted_dates)
}

# Convert calendar dates to decimal years, which is a common BEAST/TempEst format.
# Example: mid-2023 becomes approximately 2023.5.
date_to_decimal_year <- function(date) {
  year_start <- lubridate::ymd(paste0(lubridate::year(date), "-01-01"))
  next_year_start <- lubridate::ymd(paste0(lubridate::year(date) + 1, "-01-01"))
  
  lubridate::year(date) +
    as.numeric(date - year_start) / as.numeric(next_year_start - year_start)
}

# ==============================================================================
# 1. DATA LOADING
# ==============================================================================

cat("\n========== 1. DATA LOADING ==========\n\n")

###################
# 1. Load inputs
###################

# Input metadata and FASTA alignment files
metadata_file <- "data/sequences/DENV_genomic_analyses_simulated_dataset_2-3.csv"
alignment_file <- "data/sequences/DENV_genomic_analyses_simulated_dataset_2-3.fas"

check_file_exists(metadata_file)
check_file_exists(alignment_file)

metadata_all <- read.csv(metadata_file, header = TRUE, stringsAsFactors = FALSE)     # Load the complete metadata table (for Rt estimation)

# Check required metadata columns.
check_required_columns(
  metadata_all,
  required_cols = c("ID", "collection_date", "location"),
  data_name = "metadata_all"
)

# Clean important character columns.
metadata_all$ID <- trimws(metadata_all$ID)
metadata_all$collection_date <- trimws(metadata_all$collection_date)
metadata_all$location <- trimws(metadata_all$location)

# Convert collection dates once for downstream analyses
metadata_all$date <- convert_dates_safely(
  metadata_all$collection_date,
  "metadata_all$collection_date"
)

# Display basic dataset statistics
cat("Total number of cases:", nrow(metadata_all), "\n")
cat("  - H- metadata records:", sum(grepl("^H-", metadata_all$ID)), "\n")
cat("  - M- cases without genomic sequence:", sum(grepl("^M-", metadata_all$ID)), "\n")

# Separate H- and M- cases because genomic sequences can only originate
# from H- cases, whereas M- cases are never associated with sequencing data
metadata_H <- metadata_all[grepl("^H-", metadata_all$ID), ]
metadata_M <- metadata_all[grepl("^M-", metadata_all$ID), ]

cat("\nMetadata for H- cases:", nrow(metadata_H), "entries\n")
cat("Metadata for M- cases:", nrow(metadata_M), "entries\n")

# ------------------------------------------------------------------------------
# FASTA alignment loading and metadata matching
# ------------------------------------------------------------------------------
# The FASTA alignment contains the viral genome sequences used for the
# phylogenetic and phylogeographic analyses. Each FASTA header must match
# an H- case ID in the metadata table so that sampling dates and locations
# can be associated with each genome sequence.
# ------------------------------------------------------------------------------

alignment <- ape::read.dna(alignment_file, format = "fasta")

if (is.null(alignment) || length(alignment) == 0) {
  stop("FASTA file could not be read. Check that it is a valid FASTA file.")
}

# Clean sequence names.
rownames(alignment) <- trimws(rownames(alignment))

cat("\nFASTA alignment loaded:", nrow(alignment), "sequences,",
    ncol(alignment), "nucleotides\n")

# Match FASTA sequence IDs with their metadata.
fasta_ids <- rownames(alignment)
metadata_seq <- metadata_H[metadata_H$ID %in% fasta_ids, ]

# Reorder metadata to follow the exact same order as the FASTA sequences.
metadata_seq <- metadata_seq[match(fasta_ids, metadata_seq$ID), ]

if (any(is.na(metadata_seq$ID))) {
  missing_ids <- fasta_ids[is.na(metadata_seq$ID)]
  stop("Some FASTA sequences have no matching metadata: ",
       paste(missing_ids, collapse = ", "))
}

# Check metadata/alignment consistency
if (nrow(metadata_seq) != nrow(alignment)) {
  stop("The number of matched metadata rows does not match the number of FASTA sequences.")
}

cat("Sequences with matching metadata:", nrow(metadata_seq), "\n")

# Enure sequence metadata dates are in `date format
metadata_seq$date <- convert_dates_safely(
  metadata_seq$collection_date,
  "metadata_seq$collection_date"
)

# ==============================================================================
# 2. PRELIMINARY VISUALISATION OF RAW DATA
# ==============================================================================
# Summarise spatial and temporal sampling patterns before phylogenetic analyses
# ==============================================================================

cat("\n========== 2. PRELIMINARY VISUALISATION ==========\n\n")

# ------------------------------------------------------------------------------
# 2.1 Load shapefiles and prepare a regional basemap (SE Asia)
# ------------------------------------------------------------------------------

# Disable S2 to avoid topology issues when plotting this shapefile
sf::sf_use_s2(FALSE)

world_shp <- sf::st_read(
  "data/shapefiles/World_countries_shapefile.shp",
  quiet = TRUE,
  check_ring_dir = FALSE
)

# Remove the CRS for this decriptive plotting step
sf::st_crs(world_shp) <- NA

# Detect the ISO3 column name automatically: different shapefile providers
# use slightly different conventions (ISO3, ISO_A3, ADM0_A3, ...).
iso3_candidate_cols <- c("ISO3", "ISO_A3", "ADM0_A3", "iso_a3", "ADMIN_A3")
iso3_col <- iso3_candidate_cols[iso3_candidate_cols %in% names(world_shp)][1]
if (is.na(iso3_col)) {
  stop("No ISO3 column found in World_countries_shapefile.shp. ",
       "Available columns: ", paste(names(world_shp), collapse = ", "))
}
cat("ISO3 column detected in shapefile:", iso3_col, "\n")
world_shp$ISO3 <- world_shp[[iso3_col]]

# Five focal countries in this DENV dataset
focal_iso3 <- c("VNM", "THA", "KHM", "LAO", "CHN")
focal_shp  <- world_shp[world_shp$ISO3 %in% focal_iso3, ]

# Optional thin lines for international borders
borders_shp <- tryCatch(
  {
    b <- sf::st_read(
      "data/shapefiles/Only_international_borders.shp",
      quiet = TRUE,
      check_ring_dir = FALSE
    )
    
    # Remove CRS for display only
    sf::st_crs(b) <- NA
    
    b
  },
  error = function(e) NULL
)

# Cropping window for the SE Asia region (longitude/latitude coordinates)
seasia_bbox <- c(xmin = 92, ymin = 5, xmax = 125, ymax = 30)

# Country colour used across figure panels
country_colors <- c(
  VNM = "#B65A5C",  # muted red
  THA = "#5D7FA3",  # muted blue
  KHM = "#6EA87A",  # muted green
  LAO = "#B8A36A",  # muted ochre
  CHN = "#8172B3"   # muted purple
)

# ----------------------------------------------------------------------------
# 2.2 Pre-phylogeographic assessment of genomic sampling
#     Multipanel figure inspired by scientific article layouts
# ------------------------------------------------------------------------------
# Figure panels : A map, B sampling timing, C cases vss genomes, D cumulative genomes
# ------------------------------------------------------------------------------

# Per-country summary: first case, first sequence, totals, sequencing ratio
country_summary <- do.call(rbind, lapply(focal_iso3, function(c) {
  cases_c <- metadata_all[metadata_all$location == c, ]
  seqs_c  <- metadata_seq[metadata_seq$location == c, ]
  data.frame(
    location         = c,
    first_case       = if (nrow(cases_c) > 0) min(cases_c$date) else as.Date(NA),
    first_sequence   = if (nrow(seqs_c)  > 0) min(seqs_c$date)  else as.Date(NA),
    last_case        = if (nrow(cases_c) > 0) max(cases_c$date) else as.Date(NA),
    last_sequence    = if (nrow(seqs_c)  > 0) max(seqs_c$date)  else as.Date(NA),
    n_cases          = nrow(cases_c),
    n_sequences      = nrow(seqs_c),
    stringsAsFactors = FALSE
  )
}))
country_summary$lag_days <- as.numeric(country_summary$first_sequence -
                                         country_summary$first_case)
country_summary$seq_per_case <- 100 * country_summary$n_sequences /
  country_summary$n_cases
country_summary <- country_summary[order(country_summary$first_case), ]

# Compare reported cases and sequenced genomes using linear and rank-based correlations
pearson_test <- cor.test(
  country_summary$n_cases,
  country_summary$n_sequences,
  method = "pearson"
)
spearman_test <- suppressWarnings(
  cor.test(
    country_summary$n_cases,
    country_summary$n_sequences,
    method = "spearman"
  )
)
case_sequence_lm <- lm(n_sequences ~ n_cases, data = country_summary)
case_sequence_r2 <- summary(case_sequence_lm)$r.squared
overall_ratio <- 100 * sum(country_summary$n_sequences) /
  sum(country_summary$n_cases)

cat("\nCorrelation cases vs sequences (n = ", nrow(country_summary), "):\n",
    sep = "")
cat(sprintf("  Linear model R2 = %.3f\n", case_sequence_r2))
cat(sprintf("  Pearson  r   = %.3f, p = %.4f\n",
            pearson_test$estimate, pearson_test$p.value))
cat(sprintf("  Spearman rho = %.3f, p = %.4f\n",
            spearman_test$estimate, spearman_test$p.value))

cat("\nCountry summary table:\n")
print(country_summary)
write.csv(country_summary, "outputs/country_summary.csv",
          row.names = FALSE, quote = FALSE)

# Prepare map data by joining epidemiological and genomic counts to the focal country polygons.
counts_per_country <- aggregate(ID ~ location, data = metadata_all, FUN = length)
names(counts_per_country)[2] <- "n_cases"
seq_per_country <- aggregate(ID ~ location, data = metadata_seq, FUN = length)
names(seq_per_country)[2] <- "n_sequences"

map_df <- merge(focal_shp, counts_per_country,
                by.x = "ISO3", by.y = "location", all.x = TRUE)
map_df <- merge(map_df, seq_per_country,
                by.x = "ISO3", by.y = "location", all.x = TRUE)
map_df$n_cases[is.na(map_df$n_cases)] <- 0
map_df$n_sequences[is.na(map_df$n_sequences)] <- 0

# Manual SE-Asia subsetting. This uses only bounding boxes and therefore avoids
# geometry operations that can trigger S2/topology errors on country shapefiles.
bbox_overlaps_seasia <- function(geom, win) {
  bb <- sf::st_bbox(geom)
  !(bb["xmax"] < win["xmin"] || bb["xmin"] > win["xmax"] ||
      bb["ymax"] < win["ymin"] || bb["ymin"] > win["ymax"])
}
keep_world <- vapply(sf::st_geometry(world_shp),
                     bbox_overlaps_seasia, logical(1), win = seasia_bbox)
seasia_world <- world_shp[keep_world, ]

borders_crop <- NULL
if (!is.null(borders_shp)) {
  keep_borders <- vapply(sf::st_geometry(borders_shp),
                         bbox_overlaps_seasia, logical(1), win = seasia_bbox)
  borders_crop <- borders_shp[keep_borders, ]
}

# Approximate label/circle positions. These are intentionally cartographic
# positions rather than geometric centroids, to keep the map readable.
map_points <- data.frame(
  location = c("VNM", "THA", "KHM", "LAO", "CHN"),
  x = c(108.3, 101.0, 104.8, 102.8, 105.2),
  y = c(14.7, 15.8, 12.7, 18.2, 24.0),
  stringsAsFactors = FALSE
)
map_df <- merge(map_df, map_points, by.x = "ISO3", by.y = "location", all.x = TRUE)

# Circle areas are proportional to counts within each data layer. The maximum
# radii only control plotting size; the legends use the same formulas, so legend
# values and map symbols are quantitatively matched.
case_radius_max <- 0.82
seq_radius_max  <- 0.32

case_radius <- sqrt(map_df$n_cases / max(map_df$n_cases)) * case_radius_max
seq_radius  <- sqrt(map_df$n_sequences / max(map_df$n_sequences)) * seq_radius_max

# Helper function for drawing circles in map coordinates.
draw_map_circles <- function(x, y, radius, border, col = NA, lwd = 0.8) {
  symbols(
    x = x,
    y = y,
    circles = radius,
    inches = FALSE,
    add = TRUE,
    fg = border,
    bg = col,
    lwd = lwd
  )
}

# Dates used in the temporal panels.
all_dates <- seq(min(metadata_all$date), max(metadata_all$date), by = "1 day")

# Main four-panel ssampling figure
pdf("outputs/Fig1_prephylogeographic_sampling_assessment.pdf",
    width = 7.2, height = 6.65, bg = "white", useDingbats = FALSE)

layout(
  matrix(c(
    1, 2,
    1, 3,
    4, 4
  ), nrow = 3, byrow = TRUE),
  widths  = c(1.04, 1.16),
  heights = c(0.66, 1.04, 0.72)
)

par(
  oma = c(0.45, 0.25, 0.25, 0.15),
  family = "sans",
  col = "gray20",
  col.axis = "gray20",
  fg = "gray20",
  mgp = c(1.32, 0.38, 0),
  tcl = -0.18,
  cex = 0.94,
  cex.axis = 0.78,
  cex.lab = 0.82
)

draw_legend_circle <- function(x, y, r, col = NA, border = "gray20", lwd = 0.8) {
  theta <- seq(0, 2 * pi, length.out = 300)
  polygon(
    x + r * cos(theta),
    y + r * sin(theta),
    col = col,
    border = border,
    lwd = lwd,
    xpd = NA
  )
}

# ----------------------------------------------------------------------------
# Panel A
# ----------------------------------------------------------------------------

par(mar = c(0.35, 0.35, 1.05, 0.25))

plot(
  sf::st_geometry(seasia_world),
  col = "gray97",
  border = "gray75",
  lwd = 0.30,
  xlim = c(97.2, 111.0),
  ylim = c(10.2, 25.1),
  axes = FALSE,
  xlab = "",
  ylab = ""
)

if (!is.null(borders_crop)) {
  plot(sf::st_geometry(borders_crop), col = "gray65", lwd = 0.22, add = TRUE)
}

box(col = "gray25", lwd = 0.65)

for (i in seq_len(nrow(map_df))) {
  iso <- map_df$ISO3[i]
  
  outer_r <- case_radius[i]
  inner_r <- seq_radius[i]
  
  draw_map_circles(
    map_df$x[i], map_df$y[i], outer_r,
    border = adjustcolor(country_colors[iso], alpha.f = 0.85),
    col = adjustcolor(country_colors[iso], alpha.f = 0.16),
    lwd = 0.9
  )
  
  draw_map_circles(
    map_df$x[i], map_df$y[i], inner_r,
    border = adjustcolor(country_colors[iso], alpha.f = 0.98),
    col = adjustcolor(country_colors[iso], alpha.f = 0.62),
    lwd = 0.75
  )
  
  text(map_df$x[i], map_df$y[i] + outer_r + 0.55,
       labels = iso, cex = 0.64, font = 1, col = "gray25")
}

mtext("A", side = 3, line = 0.05, adj = 0, font = 2, cex = 1.18)

scale_km <- 500
scale_y <- 8.5
scale_x2 <- 110.7
scale_x1 <- scale_x2 - scale_km / (111.32 * cos(scale_y * pi / 180))

segments(scale_x1, scale_y, scale_x2, scale_y, lwd = 1.1, col = "gray15")
segments(scale_x1, scale_y - 0.08, scale_x1, scale_y + 0.08, lwd = 0.8, col = "gray15")
segments(scale_x2, scale_y - 0.08, scale_x2, scale_y + 0.08, lwd = 0.8, col = "gray15")
text(mean(c(scale_x1, scale_x2)), scale_y + 0.36, "500 km", cex = 0.56, col = "gray15")

# Reported cases legend, placed inside the map at the lower-left edge.
case_title_x <- 97.22
case_legend_x <- 97.85
case_legend_bottom <- 8
case_legend_values <- c("20,000" = 20000, "10,000" = 10000, "5,000" = 5000)
case_legend_r <- sqrt(case_legend_values / max(map_df$n_cases)) * case_radius_max
case_legend_y <- case_legend_bottom + case_legend_r
case_label_y <- c("20,000" = 8.95, "10,000" = 8.60, "5,000" = 8.32)

text(case_title_x, 10.33,
     "Reported
cases", adj = c(0, 0.5), cex = 0.54,
     col = "gray10", xpd = NA)

for (nm in names(case_legend_r)) {
  r <- case_legend_r[nm]
  y_center <- case_legend_y[nm]
  y_label <- case_label_y[nm]
  x_circle_edge <- case_legend_x +
    sqrt(pmax(r^2 - (y_label - y_center)^2, 0))
  
  draw_legend_circle(
    x = case_legend_x,
    y = y_center,
    r = r,
    col = adjustcolor("white", alpha.f = 0.80),
    border = "gray20",
    lwd = 0.75
  )
  
  segments(
    x0 = x_circle_edge,
    x1 = case_legend_x + 1.08,
    y0 = y_label,
    y1 = y_label,
    col = "gray25",
    lwd = 0.55,
    xpd = NA
  )
  
  text(
    case_legend_x + 1.20,
    y_label,
    labels = nm,
    adj = 0,
    cex = 0.54,
    col = "gray10",
    xpd = NA
  )
}

# Sequenced genomes legend, matched to the map's filled genome circles.
seq_legend_x <- 100.68
seq_legend_y <- c("100" = 9.15, "50" = 8.60, "20" = 8.20)
seq_legend_values <- c("100" = 100, "50" = 50, "20" = 20)
seq_legend_r <- sqrt(seq_legend_values / max(map_df$n_sequences)) * seq_radius_max

text(seq_legend_x - 0.45, 10.0,
     "Sequenced
genomes", adj = c(0, 0.5), cex = 0.54,
     col = "gray10", xpd = NA)

for (nm in names(seq_legend_r)) {
  draw_legend_circle(
    x = seq_legend_x,
    y = seq_legend_y[nm],
    r = seq_legend_r[nm],
    col = adjustcolor("gray10", alpha.f = 0.88),
    border = "gray10",
    lwd = 0.65
  )
  text(
    seq_legend_x + 0.34,
    seq_legend_y[nm],
    labels = nm,
    adj = 0,
    cex = 0.54,
    col = "gray10",
    xpd = NA
  )
}

# ----------------------------------------------------------------------------
# Panel B
# ----------------------------------------------------------------------------

par(mar = c(2.20, 2.20, 1.05, 0.55))

country_order_B <- c("VNM", "KHM", "CHN", "THA", "LAO")
country_summary_B <- country_summary[match(country_order_B, country_summary$location), ]

plot(
  c(min(metadata_all$date), max(metadata_all$date)),
  c(0.5, nrow(country_summary_B) + 0.5),
  type = "n",
  axes = FALSE,
  xlab = "Date",
  ylab = ""
)

for (i in seq_len(nrow(country_summary_B))) {
  ctry <- country_summary_B$location[i]
  seqs_c <- metadata_seq[metadata_seq$location == ctry, ]
  cases_c <- metadata_all[metadata_all$location == ctry, ]
  
  segments(min(cases_c$date), i, max(cases_c$date), i,
           col = "gray78", lwd = 0.8)
  
  points(seqs_c$date, rep(i, nrow(seqs_c)),
         pch = 21,
         bg = adjustcolor(country_colors[ctry], alpha.f = 0.75),
         col = adjustcolor(country_colors[ctry], alpha.f = 0.95),
         cex = 0.40)
  
  points(country_summary_B$first_case[i], i,
         pch = 21, bg = "white", col = "gray20",
         cex = 0.72, lwd = 0.75)
  
  points(country_summary_B$first_sequence[i], i,
         pch = 24, bg = country_colors[ctry], col = "gray20",
         cex = 0.74, lwd = 0.60)
}

axis(2, at = seq_len(nrow(country_summary_B)),
     labels = country_summary_B$location, las = 1, cex.axis = 0.76)

axis.Date(1,
          at = seq(min(metadata_all$date), max(metadata_all$date), by = "1 month"),
          format = "%b", cex.axis = 0.76)

box(col = "gray25", lwd = 0.65)
mtext("B", side = 3, line = 0.05, adj = 0, font = 2, cex = 1.18)

legend(
  "topleft",
  legend = c("First case", "First sequence", "Sequences"),
  pch = c(21, 24, 21),
  pt.bg = c("white", "gray60", "gray20"),
  col = "gray20",
  bty = "n",
  ncol = 1,
  cex = 0.60,
  pt.cex = 0.74,
  y.intersp = 0.72,
  x.intersp = 0.70,
  inset = c(0.01, 0.02)
)

# ----------------------------------------------------------------------------
# Panel C
# ----------------------------------------------------------------------------

par(mar = c(2.60, 3.20, 1.05, 0.45))

point_cex <- sqrt(country_summary$n_sequences) /
  sqrt(max(country_summary$n_sequences)) * 1.70

plot(
  country_summary$n_cases,
  country_summary$n_sequences,
  type = "n",
  axes = FALSE,
  xlab = "Reported cases",
  ylab = "Sequenced genomes",
  xlim = c(0, max(country_summary$n_cases) * 1.10),
  ylim = c(0, max(country_summary$n_sequences) * 1.20),
  cex.lab = 0.86
)

abline(case_sequence_lm, col = "gray45", lwd = 1.0)

points(
  country_summary$n_cases,
  country_summary$n_sequences,
  pch = 21,
  bg = adjustcolor(country_colors[country_summary$location], alpha.f = 0.75),
  col = "gray25",
  cex = point_cex,
  lwd = 0.7
)

axis(1, cex.axis = 0.68)
axis(2, las = 1, cex.axis = 0.68)
box(col = "gray25", lwd = 0.65)
mtext("C", side = 3, line = 0.05, adj = 0, font = 2, cex = 1.18)

legend(
  "bottomright",
  legend = focal_iso3,
  pch = 21,
  pt.bg = country_colors[focal_iso3],
  col = "gray25",
  bty = "n",
  ncol = 2,
  cex = 0.54,
  pt.cex = 0.66,
  y.intersp = 0.75,
  x.intersp = 0.70,
  inset = c(0.01, 0.02)
)

text(
  x = max(country_summary$n_cases) * 0.04,
  y = max(country_summary$n_sequences) * 1.10,
  labels = bquote(italic(R)^2 * " = " *
                    .(sprintf("%.2f", case_sequence_r2))),
  adj = c(0, 1),
  cex = 0.54,
  col = "gray10"
)

text(
  x = max(country_summary$n_cases) * 0.04,
  y = max(country_summary$n_sequences) * 1.00,
  labels = bquote("Pearson " * italic(r) * " = " *
                    .(sprintf("%.2f", pearson_test$estimate)) *
                    ", " * italic(P) * " = " *
                    .(sprintf("%.3f", pearson_test$p.value))),
  adj = c(0, 1),
  cex = 0.54,
  col = "gray10"
)

text(
  x = max(country_summary$n_cases) * 0.04,
  y = max(country_summary$n_sequences) * 0.90,
  labels = bquote("Spearman " * rho * " = " *
                    .(sprintf("%.2f", spearman_test$estimate)) *
                    ", " * italic(P) * " = " *
                    .(sprintf("%.3f", spearman_test$p.value))),
  adj = c(0, 1),
  cex = 0.54,
  col = "gray10"
)

# ----------------------------------------------------------------------------
# Panel D
# ----------------------------------------------------------------------------

par(mar = c(2.60, 3.20, 0.95, 0.55))

e_start_date <- as.Date("2023-06-01")
e_dates <- all_dates[all_dates >= e_start_date]

plot(
  e_dates,
  rep(0, length(e_dates)),
  type = "n",
  axes = FALSE,
  xlab = "Sampling date",
  ylab = "Cumulative sequenced\ngenomes",
  ylim = c(0, max(country_summary$n_sequences) * 1.18),
  cex.lab = 0.86
)

for (ctry in focal_iso3) {
  seq_dates_c <- metadata_seq$date[metadata_seq$location == ctry]
  cum_seq_c <- cumsum(table(factor(seq_dates_c, levels = as.character(all_dates))))
  cum_seq_e <- cum_seq_c[all_dates >= e_start_date]
  
  lines(e_dates, cum_seq_e,
        col = country_colors[ctry], lwd = 1.35)
  
  points(max(seq_dates_c), max(cum_seq_c),
         pch = 21, bg = country_colors[ctry],
         col = "gray25", cex = 0.75)
}

axis.Date(1,
          at = seq(e_start_date, max(all_dates), by = "1 month"),
          format = "%b", cex.axis = 0.68)

axis(2, las = 1, cex.axis = 0.68)
box(col = "gray25", lwd = 0.65)
mtext("D", side = 3, line = 0.05, adj = 0, font = 2, cex = 1.18)

legend(
  "top",
  legend = focal_iso3,
  col = country_colors[focal_iso3],
  lwd = 1.15,
  bty = "n",
  ncol = 5,
  cex = 0.70
)

dev.off()

cat("  -> outputs/Fig1_prephylogeographic_sampling_assessment.pdf\n")

cat("\nSampling bias summary:\n")
cat(sprintf("  Overall sequencing effort: %.2f %%\n", overall_ratio))
cat(sprintf("  Linear model R2 = %.3f\n", case_sequence_r2))
cat(sprintf("  Pearson  r   = %.3f, p = %.4f\n",
            pearson_test$estimate, pearson_test$p.value))
cat(sprintf("  Spearman rho = %.3f, p = %.4f\n",
            spearman_test$estimate, spearman_test$p.value))

# ==============================================================================
# 3. PRELIMINARY Rt ESTIMATION FROM CASE DATA
# ==============================================================================

cat("\n========== 3. Rt ESTIMATION ==========\n\n")

# Rt is estimated from all reported cases, including records without genome sequences

# Convert collection dates to epidemic days, starting at day 1
days_from_start <- lubridate::interval(
  min(metadata_all$date), 
  metadata_all$date
) %/% lubridate::days(1) + 1

total_days <- lubridate::interval(
  min(metadata_all$date), 
  max(metadata_all$date)
) %/% lubridate::days(1) + 1

# Build the daily incidence time series
daily_cases <- tabulate(days_from_start, nbins = total_days)

cat("Total number of days:", total_days, "\n")
cat("Total number of cases in daily incidence:", sum(daily_cases), "\n")
cat("Maximum daily number of cases:", max(daily_cases), "\n")

# Print basic checks for the daily incidence curve.
if (sum(daily_cases) != nrow(metadata_all)) {
  stop("Daily incidence does not sum to the total number of metadata records.")
}

# Check that daily incidence reconstruction did not losse records
top_daily_counts <- sort(table(metadata_all$collection_date), decreasing = TRUE)[1:10]

# Export the largest daily counts for reporting/QC
write.csv(
  data.frame(
    date = names(top_daily_counts),
    cases = as.integer(top_daily_counts)
  ),
  "outputs/top_daily_case_counts.csv",
  row.names = FALSE,
  quote = FALSE
)

cat("\nTop 10 daily case counts:\n")
print(top_daily_counts)

# Rt is estimated repeatedly across plausible dengue serial interval values
# Because the exact serial interval is uncertain, we repeat the Rt estimation
# several times using random values sampled from plausible ranges.

n_iterations <- 1000

# Plausible range for the mean serial interval, in days.
mean_range <- c(13, 17)

# Plausible range for the standard deviation of the serial interval, in days.
sd_range <- c(4, 6)

# Seven-day sliding windows for EpiEstim.
if (length(daily_cases) < 8) {
  stop("Not enough daily incidence data to estimate Rt with a 7-day window.")
}

# Define the start of each Rt estimation window.
# The first window starts on day 2 because EpiEstim needs previous incidence
# information to estimate transmission.
t_start <- seq(2, length(daily_cases) - 6)

# Define the end of each Rt estimation window.
# Each window lasts 7 days: day 2 to day 8, day 3 to day 9, etc.
t_end <- seq(8, length(daily_cases))

# Create an empty matrix to store Rt estimates from all iterations.
# Rows correspond to time windows, columns correspond to serial interval draws
all_Rt <- matrix(NA, nrow = length(t_start), ncol = n_iterations)

cat("\nEstimating Rt over", n_iterations, "iterations...\n")

# Fixed seed for reproductibility 
set.seed(42)

# Repeat Rt estimation under different plausible serial interval assumptions.
for (i in seq_len(n_iterations)) {
  mean_si_i <- runif(1, mean_range[1], mean_range[2]) # Randomly sample one mean serial interval from the predefined range.
  sd_si_i <- runif(1, sd_range[1], sd_range[2])       # Randomly sample one standard deviation for the serial interval.
  
  res_i <- tryCatch(                                  # Run EpiEstim.
    {                                                 # tryCatch prevents the whole script from stopping if one iteration fails.
      EpiEstim::estimate_R(
        incid = daily_cases,
        method = "parametric_si",
        config = EpiEstim::make_config(list(
          mean_si = mean_si_i,
          std_si = sd_si_i,
          t_start = t_start,
          t_end = t_end
        ))
      )
    },
    error = function(e) {
      message("Rt iteration ", i, " failed: ", e$message)
      return(NULL)
    }
  )
  
  # If the iteration succeeded, store the mean Rt estimate for each time window.  
  if (!is.null(res_i)) {
    all_Rt[, i] <- res_i$R$`Mean(R)`
  }
}

# Stop if all Rt iterations failed.
if (all(is.na(all_Rt))) {
  stop("All Rt estimation iterations failed. Check incidence data and serial interval settings.")
}

# Summarise Rt uncertainty across serial interval draws
R_median <- apply(all_Rt, 1, median, na.rm = TRUE)
R_lower <- apply(all_Rt, 1, quantile, probs = 0.025, na.rm = TRUE)
R_upper <- apply(all_Rt, 1, quantile, probs = 0.975, na.rm = TRUE)

# Assign one date to each Rt estimate.
# Dates correspond to the midpoint of each Rt estimation window.

Rt_days <- (t_start + t_end) / 2
R_dates <- min(metadata_all$date) + (Rt_days - 1)

# Store Rt estimates in a clean table for export.
Rt_results <- data.frame(
  date = R_dates,
  t_start = t_start,
  t_end = t_end,
  R_median = R_median,
  R_lower_95 = R_lower,
  R_upper_95 = R_upper
)

# Save Rt estimates for reporting and further analyses.
write.csv(
  Rt_results,
  "outputs/Rt_estimates.csv",
  row.names = FALSE,
  quote = FALSE
)

# Print the first and last valid Rt values as a quick summary.
valid_R_median <- !is.na(R_median)

if (any(valid_R_median)) {
  cat("Rt median at first estimable time point:",
      round(R_median[which(valid_R_median)[1]], 2), "\n")
  cat("Rt median at last time point:",
      round(tail(R_median[valid_R_median], 1), 2), "\n")
} else {
  warning("No valid Rt median values were produced.")
}

# ------------------------------------------------------------------------------
# Figure 4: effective reproduction number over time
# ------------------------------------------------------------------------------
# This figure shows the estimated Rt through time.
# The solid line represents the median Rt estimate.
# The shaded area represents the approximate 95% uncertainty interval.
# The red dashed horizontal line marks Rt = 1, the epidemic threshold.
# ------------------------------------------------------------------------------

pdf("outputs/Fig4_Rt_estimation.pdf", width = 12, height = 5)
par(
  oma = c(0, 0, 0, 0), mar = c(3, 4, 2, 1),
  lwd = 0.3, col = "gray30", col.axis = "gray30", fg = "gray30"
)

# Identify Rt estimates that can be plotted.
valid_rt <- !is.na(R_median) & !is.na(R_lower) & !is.na(R_upper)

# Create an empty plotting area.
plot(
  R_dates,
  R_median,
  type = "n",
  axes = FALSE,
  xlab = NA,
  ylab = NA,
  ylim = c(0, max(R_upper[valid_rt], na.rm = TRUE) * 1.1),
  main = "Effective reproduction number over time"
)

# Build the coordinates of the uncertainty polygon.
xx <- c(R_dates[valid_rt], rev(R_dates[valid_rt]))
yy <- c(R_lower[valid_rt], rev(R_upper[valid_rt]))

# Add the 95% uncertainty interval as a shaded polygon.
polygon(
  xx,
  yy,
  col = rgb(187 / 255, 187 / 255, 187 / 255, 0.3),
  border = NA
)

# Add the median Rt estimate.
lines(R_dates[valid_rt], R_median[valid_rt], lwd = 1, col = "gray30")

# Add the epidemic threshold.
# Rt > 1 indicates increasing transmission.
# Rt < 1 indicates decreasing transmission.
abline(h = 1, lty = 2, lwd = 0.5, col = "red")

# Define x-axis tick marks every two weeks.
figure4_ticks <- seq(min(R_dates), max(R_dates), by = "2 weeks")

axis(
  side = 1,
  at = figure4_ticks,
  labels = format(figure4_ticks, "%d/%m"),
  cex.axis = 0.7,
  las = 2
)

axis(side = 2, cex.axis = 0.7, las = 1)

# Add y-axis label manually for better control of spacing.
mtext("Effective reproduction number (Rt)", side = 2, line = 2.5, cex = 0.9)

dev.off()
cat("  -> outputs/Fig4_Rt_estimation.pdf\n")

# ==============================================================================
# 4. RECOMBINATION SCREENING
# ==============================================================================

cat("\n========== 4. RECOMBINATION SCREENING ==========\n\n")

# The Phi test is used to detect evidence of recombination in an alignment.
#
# Recombination can be problematic for phylogenetic analyses because standard
# tree-based methods assume that all sites in the alignment share the same
# evolutionary history.
#
# If recombination is present, different genome regions may have different
# evolutionary histories, which can make a single phylogenetic tree misleading.

# The Phi test was run externally because it was not available in the installed phangorn version

phangorn_phi_functions <- grep(
  "phi",
  ls("package:phangorn"),
  value = TRUE,
  ignore.case = TRUE
)

cat("Functions containing 'phi' in phangorn:\n")
print(phangorn_phi_functions)

# Interpretation rule:
#   - p > 0.05  : no statistically significant evidence of recombination
#   - p <= 0.05 : statistically significant evidence of recombination
#
# In this project, the p-value below was obtained externally using SplitsTree.
phi_p_value <- 0.2098 

# Store the recombination screening result in a structured table.
recombination_results <- data.frame(
  method = "Phi test",
  software = "External software required: SplitsTree or PhiPack",
  p_value = phi_p_value,
  interpretation = NA_character_
)

# Interpret the Phi test result
if (is.na(phi_p_value)) {
  
  # Case 1: the test was not run or the result was not entered.
  recombination_results$interpretation <- paste(
    "Phi test not run in R.",
    "No Phi-test function was available in the installed phangorn package.",
    "The p-value should be obtained externally and entered in the script."
  )
  
  cat("\nPhi test p-value has not been entered yet.\n")
  cat("Run the Phi test externally and replace phi_p_value <- NA_real_ with the obtained p-value.\n")
  
} else if (phi_p_value > 0.05) {
  
  # Case 2: no statistically significant recombination signal was detected
  recombination_results$interpretation <- paste(
    "No statistically significant recombination signal was detected",
    "at the 5% significance level."
  )
  
  cat("\nPhi test p-value =", round(phi_p_value, 4), "\n")
  cat("No statistically significant recombination signal was detected.\n")
} else {
  
  # Case 3: statistically significant recombination signal.
  # In this case, tree-based phylogenetic results should be interpreted carefully.
  recombination_results$interpretation <- paste(
    "A statistically significant recombination signal was detected",
    "at the 5% significance level."
  )
  
  cat("\nPhi test p-value =", round(phi_p_value, 4), "\n")
  cat("A statistically significant recombination signal was detected.\n")
}

# Export the recombination screening result
write.csv(
  recombination_results,
  "outputs/recombination_phi_test_result.csv",
  row.names = FALSE,
  quote = FALSE
)

# ==============================================================================
# 5. FILES FOR IQ-TREE 3, TEMPEST AND BEAST
# ==============================================================================

cat("\n========== 5. FILES FOR IQ-TREE 3, TEMPEST AND BEAST ==========\n\n")

# This section prepares the external input files needed for:
#   1. IQ-TREE 3  -> maximum-likelihood phylogenetic tree inference
#   2. TempEst    -> root-to-tip regression and temporal signal assessment
#   3. BEAST      -> time-calibrated phylogenetic and phylogeographic analyses

# ------------------------------------------------------------------------------
# 5.1 TempEst sampling date file
# ------------------------------------------------------------------------------

# TempEst needs a two-column tab-separated file:
#   column 1 = sequence/tip name
#   column 2 = sampling date
#
# TempEst usually expects NO header line.
# The tip names must match exactly the tip names in the IQ-TREE treefile.

tempest_dates <- metadata_seq[, c("ID", "collection_date")]

write.table(
  tempest_dates,
  file = "outputs/tempest_dates.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

cat("TempEst date file written to: outputs/tempest_dates.txt\n")

# ------------------------------------------------------------------------------
# 5.2 Clean FASTA alignment for IQ-TREE 3, TempEst and BEAST
# ------------------------------------------------------------------------------

# The original FASTA file is read earlier in the script as 'alignment'.
# Here we write a clean copy into the outputs folder

# Export a cleaned FASSTA alignment with tip labels matching the metadata
# nbcol = -1 writes each sequence on one single line.
# This is convenient for external software and avoids line-wrapping issues.

ape::write.dna(
  alignment,
  file = "outputs/DENV_alignment_cleaned.fas",
  format = "fasta",
  nbcol = -1,
  colsep = ""
)

cat("Cleaned FASTA alignment written to: outputs/DENV_alignment_cleaned.fas\n")

# ------------------------------------------------------------------------------
# 5.3 BEAST / BEAUti tip-date file
# ------------------------------------------------------------------------------

# BEAST needs sampling dates to estimate a time-calibrated phylogeny.
#
# We export two date formats:
#   - collection_date: normal calendar date, for example 2023-04-15
#   - decimal_date: decimal year format, for example 2023.286
#
# In BEAUti, use ONE format consistently.
# Decimal dates are often convenient for time-calibrated analyses.

beast_tip_dates <- data.frame(
  taxon = metadata_seq$ID,
  collection_date = metadata_seq$collection_date,
  decimal_date = date_to_decimal_year(metadata_seq$date)
)

write.table(
  beast_tip_dates,
  file = "outputs/beast_tip_dates.tsv",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

cat("BEAST tip-date file written to: outputs/beast_tip_dates.tsv\n")

# ------------------------------------------------------------------------------
# 5.4 BEAST / BEAUti discrete location trait file
# ------------------------------------------------------------------------------

# For a discrete phylogeographic BEAST analysis, each sequence needs a location.
#
# Here, the location is the country ISO3 code:
#   VNM = Vietnam
#   THA = Thailand
#   KHM = Cambodia
#   LAO = Laos
#   CHN = China
#
# In BEAUti, this file can be used to define a discrete trait named "location".
# Without this trait, BEAST can estimate a dated phylogeny but not viral movement between countries.

beast_discrete_location_traits <- data.frame(
  taxon = metadata_seq$ID,
  location = metadata_seq$location
)

write.table(
  beast_discrete_location_traits,
  file = "outputs/beast_discrete_location_traits.tsv",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

cat("BEAST discrete location trait file written to: outputs/beast_discrete_location_traits.tsv\n")

# ------------------------------------------------------------------------------
# 5.5 IQ-TREE 3 command templates
# ------------------------------------------------------------------------------

# IQ-TREE 3 is run externally. the command template is written to outputs/iqtree_commands.txt
#
# Main output file used later: outputs/DENV_ML.treefile

iqtree_commands <- c(
  "# IQ-TREE 3 command template",
  "# =========================",
  "",
  "# Open a Terminal in the project root folder.",
  "# The project root is the folder containing data/ and outputs/.",
  "",
  "# Then run the recommended command below:",
  "",
  "# Recommended IQ-TREE 3 command:",
  "iqtree3 -s outputs/DENV_alignment_cleaned.fas -m MFP -B 1000 --prefix outputs/DENV_ML",
  "",
  "# Command explanation:",
  "#   iqtree3",
  "#     Calls IQ-TREE version 3.",
  "#",
  "#   -s outputs/DENV_alignment_cleaned.fas",
  "#     Input FASTA alignment.",
  "#",
  "#   -m MFP",
  "#     Runs ModelFinder Plus.",
  "#     IQ-TREE tests substitution models and uses the best model automatically.",
  "#",
  "#   -B 1000",
  "#     Runs 1000 ultrafast bootstrap replicates.",
  "#     This gives branch support values.",
  "#",
  "#   --prefix outputs/DENV_ML",
  "#     Sets the output file prefix.",
  "#     This keeps IQ-TREE output files grouped under the same prefix.",
  "",
  "# Main output files produced by IQ-TREE 3:",
  "#   outputs/DENV_ML.iqtree    = full IQ-TREE report, including selected model",
  "#   outputs/DENV_ML.treefile  = maximum-likelihood tree used for TempEst and plotting",
  "#   outputs/DENV_ML.log       = IQ-TREE log file",
  "#   outputs/DENV_ML.contree   = consensus tree with bootstrap support",
  "",
  "# Alternative step-by-step workflow:",
  "",
  "# 1. Run ModelFinder only:",
  "iqtree3 -s outputs/DENV_alignment_cleaned.fas -m MF --prefix outputs/DENV_modeltest",
  "",
  "# 2. Then run tree inference manually with the selected model:",
  "# Replace <BEST_MODEL> with the best model reported in outputs/DENV_modeltest.iqtree",
  "iqtree3 -s outputs/DENV_alignment_cleaned.fas -m <BEST_MODEL> -B 1000 --prefix outputs/DENV_ML"
)

writeLines(iqtree_commands, "outputs/iqtree_commands.txt")

cat("IQ-TREE 3 command templates written to: outputs/iqtree_commands.txt\n")
cat("Recommended IQ-TREE 3 treefile will be: outputs/DENV_ML.treefile\n")

# ==============================================================================
# 6. ML TREE VISUALISATION WITH FIGTREE
# ==============================================================================

cat("\n========== 6. ML TREE VISUALISATION WITH FIGTREE ==========\n\n")

# The ML tree is visualised externally in FigTree

treefile_path <- "outputs/DENV_ML.treefile"

if (file.exists(treefile_path)) {
  cat("IQ-TREE treefile found:", treefile_path, "\n")
} else {
  cat("IQ-TREE treefile not found yet. Run the command in outputs/iqtree_commands.txt first.\n")
}

# ==============================================================================
# 7. EXTERNAL SOFTWARE PROTOCOL
# ==============================================================================

cat("\n========== 7. EXTERNAL SOFTWARE PROTOCOL ==========\n\n")

external_protocol <- c(
  "EXTERNAL SOFTWARE STEPS FOR THE DENV GENOMIC ANALYSIS",
  "==================================================",
  "",
  "1. SplitsTree4 / Phi test",
  "   - Open outputs/DENV_alignment_cleaned.fas in SplitsTree4.",
  "   - Run the Phi test for recombination.",
  "   - Report the p-value.",
  "   - If p > 0.05, there is no statistically significant evidence of recombination.",
  "   - If p <= 0.05, recombination may violate the assumption of a single tree-like history.",
  "",
  "2. IQ-TREE",
  "   - Open a terminal in the project root folder.",
  "   - Run the commands written in outputs/iqtree_commands.txt.",
  "   - Use the best-fit substitution model selected by ModelFinder.",
  "   - The ML tree file produced by the commands is outputs/DENV_ML.treefile.",
  "",
  "3. TempEst",
  "   - Open the IQ-TREE treefile in TempEst.",
  "   - Import sampling dates from outputs/tempest_dates.txt.",
  "   - Check root-to-tip divergence against sampling time.",
  "   - A positive root-to-tip relationship was used as evidence of temporal signal.",
  "",
  "4. BEAUti / BEAST: time-calibrated phylogenetics",
  "   - Import outputs/DENV_alignment_cleaned.fas into BEAUti.",
  "   - Import or enter tip dates from outputs/beast_tip_dates.tsv.",
  "   - Use decimal dates or calendar dates consistently.",
  "   - Select an appropriate substitution model based on IQ-TREE/ModelFinder.",
  "   - Choose a molecular clock model after considering the TempEst result.",
  "   - Choose a tree prior appropriate for an epidemic dataset.",
  "",
  "5. BEAUti / BEAST: phylogeography",
  "   - To make the analysis phylogeographic, add a location trait for each taxon.",
  "   - Use outputs/beast_discrete_location_traits.tsv for discrete country-level phylogeography.",
  "   - The trait column is named 'location' and contains VNM, THA, KHM, LAO, or CHN.",
  "   - In BEAUti, configure a discrete trait / discrete phylogeographic model on this location trait.",
  "   - Without this trait, this analysis remains dated phylogennetics rather than discrete phylogeography.",
  "   - If the course requires continuous phylogeography, create a separate coordinate table",
  "     with taxon/longitude/latitude columns; this script currently exports only discrete country states.",
  "",
  "6. Tracer",
  "   - Open the BEAST .log file in Tracer.",
  "   - Check ESS values; ESS > 200 was used as a practical convergence threshold.",
  "   - Inspect whether posterior, likelihood, clock rate, and tree parameters mix well.",
  "",
  "7. TreeAnnotator",
  "   - Summarise the posterior tree distribution after removing an appropriate burn-in.",
  "   - Produce a maximum clade credibility tree.",
  "",
  "8. FigTree",
  "   - Open the MCC tree.",
  "   - Display node ages, posterior probabilities, and location states if available.",
  "   - Use the location trait to colour branches or tips for phylogeographic interpretation."
)

writeLines(external_protocol, "outputs/external_software_protocol.txt")
cat("External software protocol written to: outputs/external_software_protocol.txt\n")

cat("\n========== ANALYSIS SCRIPT COMPLETED ==========\n")