# ==============================================================================
# Example.R
# Reproducible Simulated Example: Planar (2D) vs. Linear Network Comparison
# ==============================================================================
# Accompanying code for the paper:
#   Mrkvicka, T., & Kraft, S. (2026). A Permutation-Based Intensity Test for
#   Comparing Spatial Point Processes: Motorcycle Crash Locations and Perceived
#   Road Risks on a Linear Network.
#
# PURPOSE OF THIS SCRIPT:
#   This self-contained script demonstrates how to use both functions in this
#   repository (Comparison_of_2pp_2D.R and Comparison_of_2pp_LN.R) without needing any external
#   data files. It simulates a realistic road-safety scenario on a 10 x 10 km
#   road network mimicking the case study in Section 5 of the paper:
#
#     - Pattern 1 (X1 / LPP1): "Objective Crashes" (n1 = 600 events)
#     - Pattern 2 (X2 / LPP2): "Perceived Road Risks" (n2 = 300 events)
#
#   We plant two localized anomalies along the road network:
#     1. UNRECOGNIZED HAZARD (North-West, around x = 2500 m, y = 7800 m):
#        A dangerous road section where actual crashes occur frequently, but
#        drivers/riders rarely perceive it as hazardous (f1 > f2 -> RED).
#     2. OVERESTIMATED RISK (South-East, around x = 7800 m, y = 2500 m):
#        A visually intimidating road section where riders frequently report
#        fear/hazard, but actual crashes are rare (f1 < f2 -> GREEN/BLUE).
#
#   We then run:
#     - Step A: The Planar 2D Test (test_permutation_kde_env from Crashes.R),
#               which identifies broad regional zones of mismatch.
#     - Step B: The Network-Constrained Test (test_permutation_kde_env_linnet
#               from CrashesLN.R), which pinpoints the exact road segments.
# ==============================================================================

# --- 1. Load required packages and source repository functions ----------------
required_pkgs <- c("spatstat", "spatstat.linnet", "GET", "foreach", "doParallel", "doRNG")
missing_pkgs  <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop(sprintf(
    "Missing required R packages: %s
Please install them first using:
install.packages(c(%s))",
    paste(missing_pkgs, collapse = ", "),
    paste(sprintf('"%s"', missing_pkgs), collapse = ", ")
  ))
}

source("Comparison_of_2pp_2D.R")
source("Comparison_of_2pp_LN.R")

set.seed(2026)

# --- 2. Construct a 10 km x 10 km road network --------------------------------
# We rescale the built-in simplenet road network from [0, 1]^2 to a realistic
# 10,000 m x 10,000 m (10 km x 10 km) regional window in meters.
L_road <- affine(simplenet, mat = diag(c(10000, 10000)))
W_area <- Window(L_road)

cat("======================================================================
")
cat("Simulating Objective Crashes vs. Perceived Risks on a 10x10 km network
")
cat("======================================================================
")

# --- 3. Define spatial intensity functions on the road network ----------------
# Baseline intensity is uniform along the road network, plus two Gaussian hotspots:
#   - Hotspot A (North-West: x = 2500, y = 7800): Excess objective crashes (f1 > f2)
#   - Hotspot B (South-East: x = 7800, y = 2500): Excess perceived risk (f1 < f2)

intensity_crashes_fn <- function(x, y, seg, tp) {
  base_rate <- 1.0
  hotspot_A <- 12.0 * exp(-((x - 2500)^2 + (y - 7800)^2) / (2 * 650^2))
  base_rate + hotspot_A
}

intensity_perceived_fn <- function(x, y, seg, tp) {
  base_rate <- 1.0
  hotspot_B <- 12.0 * exp(-((x - 7800)^2 + (y - 2500)^2) / (2 * 650^2))
  base_rate + hotspot_B
}

# Generate exact fixed sample sizes (n1 = 600 crashes, n2 = 300 perceived risks)
n1_crashes   <- 600
n2_perceived <- 300

# Create pixel images of the unnormalized intensities on the network to sample from
im_crashes   <- as.linim(intensity_crashes_fn, L = L_road)
im_perceived <- as.linim(intensity_perceived_fn, L = L_road)

# Normalize to integrate to n1 and n2 along the network, then draw exact counts
LPP_crashes   <- rpoislpp(im_crashes / integral(im_crashes) * n1_crashes * 1.5, L = L_road)
while (npoints(LPP_crashes) < n1_crashes) {
  LPP_crashes <- superimpose(LPP_crashes, rpoislpp(im_crashes / integral(im_crashes) * 200, L = L_road))
}
LPP_crashes <- LPP_crashes[sample.int(npoints(LPP_crashes), n1_crashes)]

LPP_perceived <- rpoislpp(im_perceived / integral(im_perceived) * n2_perceived * 1.5, L = L_road)
while (npoints(LPP_perceived) < n2_perceived) {
  LPP_perceived <- superimpose(LPP_perceived, rpoislpp(im_perceived / integral(im_perceived) * 100, L = L_road))
}
LPP_perceived <- LPP_perceived[sample.int(npoints(LPP_perceived), n2_perceived)]

# Extract the 2D planar point patterns (ignoring the network) for Step A
X_crashes   <- as.ppp(LPP_crashes)
X_perceived <- as.ppp(LPP_perceived)

# --- 4. Visualize the simulated datasets --------------------------------------
par(mfrow = c(1, 2), mar = c(1, 1, 3, 1))
plot(LPP_crashes,
     main = sprintf("1. Objective Crashes (n1 = %d)
[Unrecognized hotspot in NW]", n1_crashes),
     cols = "firebrick", pch = 16, cex = 0.6, use.marks = FALSE)
points(2500, 7800, pch = 1, cex = 4, lwd = 2, col = "darkred")

plot(LPP_perceived,
     main = sprintf("2. Perceived Risks (n2 = %d)
[Overestimated fear hotspot in SE]", n2_perceived),
     cols = "dodgerblue4", pch = 16, cex = 0.6, use.marks = FALSE)
points(7800, 2500, pch = 1, cex = 4, lwd = 2, col = "navy")
par(mfrow = c(1, 1))


# ==============================================================================
# STEP A: PLANAR (2D) PERMUTATION KDE TEST
# ==============================================================================
cat("
--- Step A: Running Planar (2D) Permutation KDE Test ---
")

# Bandwidth h = 800 m on a 64x64 grid (using 499 permutations for fast demo;
# increase to 1999+ for publication results)
res_planar <- test_permutation_kde_env(
  X1 = X_crashes,
  X2 = X_perceived,
  h = 800,
  n_permutations = 499,
  fdr_level = 0.05,
  dimyx = c(64, 64)
)

# Plot planar test results (observed difference, envelope bounds, and significant regions)
plot(res_planar, what = c("obs", "lo", "hi", "lo.sign", "hi.sign"))


# ==============================================================================
# STEP B: NETWORK-CONSTRAINED PERMUTATION KDE TEST
# ==============================================================================
cat("
--- Step B: Running Network-Constrained Permutation KDE Test ---
")

# Bandwidth h = 450 m on a 128x128 network grid (using 499 parallel permutations
# for fast demo; increase to 1999+ for publication results)
res_network <- test_permutation_kde_env_linnet(
  LPP1 = LPP_crashes,
  LPP2 = LPP_perceived,
  h = 450,
  n_permutations = 499,
  fdr_level = 0.05,
  dimyx = c(128, 128),
  num_cores = 2
)

# Plot network-constrained test results
plot(res_network, what = c("obs", "lo", "hi", "lo.sign", "hi.sign"))

cat("
======================================================================
")
cat("INTERPRETATION GUIDE FOR GEOGRAPHERS:
")
cat("  * RED areas ('hi.sign', obs > hi): Objective crash density significantly
")
cat("    exceeds perceived risk (f1 > f2). These are 'hidden hazards' where
")
cat("    road users underestimate true accident risk (North-West road sections).
")
cat("  * GREEN/BLUE areas ('lo.sign', obs < lo): Perceived risk density
")
cat("    significantly exceeds actual crash records (f1 < f2). These are places
")
cat("    where road users overestimate danger (South-East road sections).
")
cat("======================================================================
")

