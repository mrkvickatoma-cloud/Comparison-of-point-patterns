# ==============================================================================
# Comparison_of_2pp_2D.R
# Planar (2D) Permutation-Based Kernel Density Difference Test with Global FDR
# ==============================================================================
# Accompanying code for the paper:
#   Mrkvicka, T., & Kraft, S. (2026). A Permutation-Based Intensity Test for
#   Comparing Spatial Point Processes: Motorcycle Crash Locations and Perceived
#   Road Risks on a Linear Network.
#
# DESCRIPTION FOR GEOGRAPHERS AND APPLIED USERS:
#   This script provides the function test_permutation_kde_env(), which
#   compares the spatial distribution (normalized intensity / probability density)
#   of two planar point patterns (e.g., actual traffic crashes vs. perceived
#   hazardous locations) across a 2D study area.
#
#   Key advantages over quadrat-based tests (SPPT, Fisher's exact, bootstrap):
#     1. Operates directly on continuous Kernel Density Estimation (KDE) surfaces,
#        avoiding arbitrary areal partitioning (MAUP).
#     2. Handles unequal sample sizes (n1 != n2) by comparing normalized densities
#        and generating the null distribution via exact random relabeling.
#     3. Controls for multiple testing across the entire map using a two-sided
#        Global False Discovery Rate (FDR) envelope (GET package).
#
# ARGUMENTS OF test_permutation_kde_env():
#   X1             : Planar point pattern (ppp object from 'spatstat') for group 1.
#   X2             : Planar point pattern (ppp object from 'spatstat') for group 2.
#   h              : Numeric smoothing bandwidth (sigma of Gaussian kernel) in
#                    map units (e.g., meters). If NULL, estimated automatically
#                    via likelihood cross-validation (bw.ppl) on pooled data.
#   n_permutations : Number of random-relabeling permutations (default = 1999).
#   fdr_level      : Nominal False Discovery Rate level alpha (default = 0.05).
#   dimyx          : Grid resolution c(ny, nx) for KDE evaluation (default = c(64, 64)).
#   verbose        : Logical; if TRUE, prints progress bar and summary.
# ==============================================================================

library(spatstat)
library(GET)

test_permutation_kde_env <- function(X1,
                                     X2,
                                     h = NULL,
                                     n_permutations = 1999,
                                     fdr_level = 0.05,
                                     dimyx = c(64, 64),
                                     verbose = TRUE) {
  
  # --- Step 0: Input validation and preprocessing ---
  if (!inherits(X1, "ppp") || !inherits(X2, "ppp")) {
    stop("Both X1 and X2 must be planar point pattern ('ppp') objects from the 'spatstat' package.")
  }
  
  # Remove any attached attribute columns (marks) so that pooling works cleanly
  X1 <- unmark(X1)
  X2 <- unmark(X2)
  
  n1 <- npoints(X1)
  n2 <- npoints(X2)
  if (n1 < 2 || n2 < 2) {
    stop("Both point patterns X1 and X2 must contain at least 2 points.")
  }
  
  # Pool all points into a single combined pattern under H0
  X_pool <- superimpose(X1, X2)
  
  # Automatic bandwidth selection via likelihood cross-validation if h is not provided
  if (is.null(h)) {
    h <- as.numeric(bw.ppl(X_pool))
    if (verbose) {
      cat(sprintf("[Planar KDE Test] Automatically selected bandwidth h = %.2f (via bw.ppl)
", h))
    }
  } else if (verbose) {
    cat(sprintf("[Planar KDE Test] Using user-specified bandwidth h = %.2f
", h))
  }
  
  # --- Step 1: Compute the observed normalized KDE difference T_0(u) ---
  # Edge correction is omitted because boundary truncation affects both
  # normalized surfaces symmetrically on the same window W under H0.
  kde1 <- density(X1, sigma = h, dimyx = dimyx, edge = FALSE)
  kde2 <- density(X2, sigma = h, dimyx = dimyx, edge = FALSE)
  
  # Normalize intensities into spatial probability densities (integrating to 1)
  T_obs <- (kde1 / n1) - (kde2 / n2)
  
  # Extract pixel grid coordinates (spatstat matrices $v are stored column-wise:
  # y-coordinates correspond to rows and vary fastest, x-coordinates to columns)
  grid_xy <- expand.grid(y = T_obs$yrow, x = T_obs$xcol)
  
  # Flatten the observed difference matrix and filter out pixels outside the study window (NA)
  v_obs <- as.vector(T_obs$v)
  noNA_id <- which(!is.na(v_obs))
  obs_valid <- v_obs[noNA_id]
  
  if (verbose) {
    cat(sprintf(
      "[Planar KDE Test] Running %d permutations across %d valid grid pixels (n1 = %d, n2 = %d)...
",
      n_permutations, length(noNA_id), n1, n2
    ))
  }
  
  # Pre-allocate matrix for permuted test statistics
  sim_m_valid <- matrix(0, nrow = length(noNA_id), ncol = n_permutations)
  
  # --- Step 2: Random-relabeling permutations under H0 ---
  if (verbose) {
    pb <- txtProgressBar(min = 0, max = n_permutations, style = 3)
  }
  
  for (k in seq_len(n_permutations)) {
    # Randomly draw n1 points without replacement for group 1; the remaining n2 go to group 2
    idx1 <- sample.int(n1 + n2, size = n1, replace = FALSE)
    
    kde1_perm <- density(X_pool[idx1], sigma = h, dimyx = dimyx, edge = FALSE)
    kde2_perm <- density(X_pool[-idx1], sigma = h, dimyx = dimyx, edge = FALSE)
    
    T_perm <- (kde1_perm / n1) - (kde2_perm / n2)
    
    # Store only valid (in-window) pixels as column k
    sim_m_valid[, k] <- as.vector(T_perm$v)[noNA_id]
    
    if (verbose && (k %% max(1, floor(n_permutations / 50)) == 0 || k == n_permutations)) {
      setTxtProgressBar(pb, k)
    }
  }
  if (verbose) close(pb)
  
  # --- Step 3: Build spatial curve_set object for the GET package ---
  r_df <- data.frame(
    x = grid_xy$x[noNA_id],
    y = grid_xy$y[noNA_id],
    width = T_obs$xstep,
    height = T_obs$ystep
  )
  
  cset_spatial <- create_curve_set(list(
    r = r_df,
    obs = obs_valid,
    sim_m = sim_m_valid
  ))
  
  # --- Step 4: Compute two-sided Global False Discovery Rate (FDR) Envelope ---
  res_fdr_spatial <- fdr_envelope(
    cset_spatial,
    alpha = fdr_level,
    alternative = "two.sided"
  )
  
  if (verbose) {
    n_high <- sum(res_fdr_spatial$obs > res_fdr_spatial$hi)
    n_low  <- sum(res_fdr_spatial$obs < res_fdr_spatial$lo)
    cat(sprintf(
      "[Planar KDE Test] Done! Significant pixels: %d above (f1 > f2), %d below (f1 < f2).
",
      n_high, n_low
    ))
  }
  
  return(res_fdr_spatial)
}


# ==============================================================================
# EMPIRICAL WORKFLOW TEMPLATE (Czech Republic Motorcycle Crashes - Planar 2D)
# ==============================================================================
# Note: The code block below reproduces the empirical planar analysis from
# Section 5.2 of the paper when you have the local CSV datasets in your working
# directory. It is wrapped in if (FALSE) so that source("Crashes.R") only
# loads the function test_permutation_kde_env() without throwing file errors.
# For a self-contained runnable example on simulated data, run Example.R.
# ==============================================================================
if (FALSE) {
  library(sf)
  library(rnaturalearth)
  # For high-resolution boundaries, install rnaturalearthhires:
  # install.packages("rnaturalearthhires", repos = "https://ropensci.r-universe.dev")
  library(rnaturalearthhires)
  
  # 1. Load point data (coordinates X, Y in S-JTSK / EPSG:5514, meters)
  crashes_df   <- read.csv("MOtonehody21_24_2.csv", header = TRUE, sep = ";", dec = ",")
  perceived_df <- read.csv("Motovyzkum_KGE.csv", header = TRUE, sep = ";", dec = ",")
  
  # 2. Obtain Czech Republic boundary and project to S-JTSK (EPSG:5514)
  cz_sf        <- ne_countries(country = "Czechia", scale = "large", returnclass = "sf")
  cz_projected <- st_transform(cz_sf, crs = 5514)
  cz_owin      <- as.owin(cz_projected)
  
  # 3. Create planar point pattern (ppp) objects clipped to the Czech boundary
  X_crashes   <- ppp(crashes_df$X, crashes_df$Y, window = cz_owin)[cz_owin]
  X_perceived <- ppp(perceived_df$X, perceived_df$Y, window = cz_owin)[cz_owin]
  
  # 4. Plot both point patterns side by side
  par(mfrow = c(1, 2), mar = c(1, 1, 3, 1))
  plot(X_crashes, main = "Objective crash locations (Police)",
       cols = "black", pch = 16, cex = 0.5, use.marks = FALSE)
  plot(X_perceived, main = "Perceived hazardous locations (Survey)",
       cols = "black", pch = 16, cex = 0.5, use.marks = FALSE)
  par(mfrow = c(1, 1))
  
  # 5. Estimate common bandwidth h via likelihood cross-validation
  h_est <- mean(c(as.numeric(bw.ppl(X_crashes)), as.numeric(bw.ppl(X_perceived))))
  
  # 6. Run the planar permutation KDE test (19,999 permutations on a 64x64 grid)
  res_2d <- test_permutation_kde_env(
    X1 = X_crashes,
    X2 = X_perceived,
    h = h_est,
    n_permutations = 19999,
    fdr_level = 0.05,
    dimyx = c(64, 64)
  )
  
  # 7. Summarize and plot the results
  cat("Total significant pixels:", sum(res_2d$obs > res_2d$hi) + sum(res_2d$obs < res_2d$lo), "
")
  plot(res_2d, what = c("obs", "lo", "hi", "lo.sign", "hi.sign"))
}
