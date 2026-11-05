# ==============================================================================
# Comparison_of_2pp_LN.R
# Linear Network Permutation-Based Kernel Density Difference Test with Global FDR
# ==============================================================================
# Accompanying code for the paper:
#   Mrkvicka, T., & Kraft, S. (2026). A Permutation-Based Intensity Test for
#   Comparing Spatial Point Processes: Motorcycle Crash Locations and Perceived
#   Road Risks on a Linear Network.
#
# DESCRIPTION FOR GEOGRAPHERS AND APPLIED USERS:
#   Many geographical events (traffic crashes, street crimes, infrastructure
#   faults) occur strictly along a linear network (roads, rivers, utility lines).
#   Applying standard 2D planar KDE spreads probability into empty off-road
#   areas and only identifies broad regional zones.
#
#   This script provides test_permutation_kde_env_linnet(), which restricts
#   both spatial point patterns (lpp objects) and kernel evaluation directly
#   to the underlying linear network L. It combines:
#     1. Fast 2D Euclidean kernel smoothing restricted to the 1D network pixels
#        (distance = "euclidean" in spatstat.linnet).
#     2. Multi-core parallel computing (foreach, doParallel, doRNG).
#     3. Two-sided Global False Discovery Rate (FDR) envelopes (GET package).
#
# ARGUMENTS OF test_permutation_kde_env_linnet():
#   LPP1           : Point pattern on a linear network (lpp object) for group 1.
#   LPP2           : Point pattern on a linear network (lpp object) for group 2.
#   h              : Numeric smoothing bandwidth (sigma of Gaussian kernel) in
#                    network coordinate units (e.g., meters).
#   n_permutations : Number of random-relabeling permutations (default = 1999).
#   fdr_level      : Nominal False Discovery Rate level alpha (default = 0.05).
#   dimyx          : Raster discretization resolution c(ny, nx) intersected with
#                    the linear network (default = c(256, 256)).
#   num_cores      : Number of CPU cores for parallel execution. If NULL, uses
#                    all available cores minus 1.
#   distance       : Kernel type on the network: "euclidean" (default, fast) or
#                    "path" (exact shortest-path kernel for small networks).
#   verbose        : Logical; if TRUE, prints status messages.
# ==============================================================================

library(spatstat)
library(spatstat.linnet)
library(GET)
library(foreach)
library(doParallel)
library(doRNG)

test_permutation_kde_env_linnet <- function(LPP1,
                                            LPP2,
                                            h,
                                            n_permutations = 1999,
                                            fdr_level = 0.05,
                                            dimyx = c(256, 256),
                                            num_cores = NULL,
                                            distance = "euclidean",
                                            verbose = TRUE) {
  
  # --- Step 0: Input validation and parallel cluster setup ---
  if (!inherits(LPP1, "lpp") || !inherits(LPP2, "lpp")) {
    stop("Both LPP1 and LPP2 must be linear network point pattern ('lpp') objects from 'spatstat.linnet'.")
  }
  if (missing(h) || is.null(h) || !is.numeric(h) || h <= 0) {
    stop("Please provide a positive numeric smoothing bandwidth 'h' (in network coordinate units).")
  }
  
  # Strip marks (attributes) to avoid conflicts when pooling patterns
  LPP1 <- unmark(LPP1)
  LPP2 <- unmark(LPP2)
  
  n1 <- npoints(LPP1)
  n2 <- npoints(LPP2)
  if (n1 < 2 || n2 < 2) {
    stop("Both point patterns LPP1 and LPP2 must contain at least 2 points.")
  }
  
  # Pool observations on the shared linear network
  LPP_pool <- superimpose(LPP1, LPP2)
  
  # Determine number of CPU cores for parallel processing
  if (is.null(num_cores)) {
    detected <- parallel::detectCores()
    num_cores <- if (is.na(detected)) 1L else max(1L, detected - 1L)
  }
  
  # --- Step 1: Compute observed normalized KDE difference T_0(u) on the network ---
  kde1 <- density(LPP1, sigma = h, distance = distance, dimyx = dimyx)
  kde2 <- density(LPP2, sigma = h, distance = distance, dimyx = dimyx)
  
  T_obs <- (kde1 / n1) - (kde2 / n2)
  
  # Extract coordinates of evaluation pixels along the discretized network
  grid_xy <- expand.grid(y = T_obs$yrow, x = T_obs$xcol)
  v_obs <- as.vector(T_obs$v)
  noNA_id <- which(!is.na(v_obs))
  obs_valid <- v_obs[noNA_id]
  
  if (verbose) {
    cat(sprintf(
      "[Network KDE Test] Discretized network into %d evaluation segments (grid: %dx%d, h = %.2f).
",
      length(noNA_id), dimyx[1], dimyx[2], h
    ))
    cat(sprintf(
      "[Network KDE Test] Running %d parallel permutations across %d CPU core(s)...
",
      n_permutations, num_cores
    ))
  }
  
  # --- Step 2: Parallel random-relabeling permutations under H0 ---
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  # Ensure the parallel cluster is always stopped cleanly, even if interrupted
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  # %dorng% guarantees reproducible, independent random number streams across workers
  sim_m_valid <- foreach(
    k = seq_len(n_permutations),
    .combine = cbind,
    .packages = c("spatstat.geom", "spatstat.linnet")
  ) %dorng% {
    
    idx1 <- sample.int(n1 + n2, size = n1, replace = FALSE)
    
    LPP1_perm <- LPP_pool[idx1]
    LPP2_perm <- LPP_pool[-idx1]
    
    kde1_perm <- density(LPP1_perm, sigma = h, distance = distance, dimyx = dimyx)
    kde2_perm <- density(LPP2_perm, sigma = h, distance = distance, dimyx = dimyx)
    
    T_perm <- (kde1_perm / n1) - (kde2_perm / n2)
    
    # Return vector of differences at valid network pixels (combined as columns)
    as.vector(T_perm$v)[noNA_id]
  }
  
  # --- Step 3: Build spatial curve_set object with network pixel geometry ---
  r_df <- data.frame(
    x = grid_xy$x[noNA_id],
    y = grid_xy$y[noNA_id],
    width = T_obs$xstep,
    height = T_obs$ystep
  )
  
  cset_spatial <- create_curve_set(list(
    r = r_df,
    obs = obs_valid,
    sim_m = as.matrix(sim_m_valid)
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
      "[Network KDE Test] Done! Significant network pixels: %d above (f1 > f2), %d below (f1 < f2).
",
      n_high, n_low
    ))
  }
  
  return(res_fdr_spatial)
}


# ==============================================================================
# EMPIRICAL WORKFLOW TEMPLATE (Czech Republic Road Network Analysis)
# ==============================================================================
# Note: The code block below reproduces the national road network analysis from
# Section 5.3 of the paper when you have the local CSV datasets in your working
# directory. It is wrapped in if (FALSE) so that source("CrashesLN.R") only
# loads the function test_permutation_kde_env_linnet() without throwing errors.
# For a self-contained runnable example on simulated data, run Example.R.
# ==============================================================================
if (FALSE) {
  # 1. Load network nodes (vertices), road edges, and point datasets (in S-JTSK meters)
  nodes_df     <- read.csv("01_krizovatky_ID.csv", header = TRUE, sep = ";")
  edges_df     <- read.csv("02_krizovtaky_spojnice.csv", header = TRUE, sep = ";")
  crashes_df   <- read.csv("MOtonehody21_24_2.csv", header = TRUE, sep = ";", dec = ",")
  perceived_df <- read.csv("Motovyzkum_KGE.csv", header = TRUE, sep = ";", dec = ",")
  
  # 2. Define bounding observation window W and build the linear network L
  W <- owin(xrange = c(-904285, -431835), yrange = c(-1227046, -935546))
  vertices <- ppp(nodes_df$XO_sjtsk, nodes_df$YO_sjtsk, window = W)
  edges    <- as.matrix(edges_df)
  L        <- linnet(vertices = vertices, edges = edges)
  
  # 3. Project both point patterns onto the linear road network (creating lpp objects)
  X_crashes     <- ppp(crashes_df$X, crashes_df$Y, window = W)
  LPP_crashes   <- lpp(X_crashes, L)
  
  X_perceived   <- ppp(perceived_df$X, perceived_df$Y, window = W)
  LPP_perceived <- lpp(X_perceived, L)
  
  # 4. Plot both network-constrained point patterns side by side
  par(mfrow = c(1, 2), mar = c(1, 1, 3, 1))
  plot(LPP_crashes, main = "Objective crash locations on network",
       cols = "black", pch = 16, cex = 0.5, use.marks = FALSE)
  plot(LPP_perceived, main = "Perceived hazardous locations on network",
       cols = "black", pch = 16, cex = 0.5, use.marks = FALSE)
  par(mfrow = c(1, 1))
  
  # 5. Run the network-constrained permutation KDE test
  # Note: For quick testing use n_permutations = 199 or 1999; for the full
  # national paper reproduction (~27,000 segments) use n_permutations = 49999.
  sigma_bw <- 1250
  res_net <- test_permutation_kde_env_linnet(
    LPP1 = LPP_crashes,
    LPP2 = LPP_perceived,
    h = sigma_bw,
    n_permutations = 1999,
    fdr_level = 0.05,
    dimyx = c(256, 256)
  )
  
  # 6. Print counts of significant network pixels and plot the envelope results
  cat("Significant network pixels (High / f1 > f2):", sum(res_net$obs > res_net$hi), "
")
  cat("Significant network pixels (Low  / f1 < f2):", sum(res_net$obs < res_net$lo), "
")
  
  plot(res_net, what = c("obs", "lo", "hi", "lo.sign", "hi.sign"))
}
