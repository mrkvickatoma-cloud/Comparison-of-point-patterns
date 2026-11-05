# Permutation-Based Intensity Test for Comparing Spatial Point Processes (Planar 2D & Linear Networks)

Accompanying R code for the paper:  
**Mrkvička, T., & Kraft, S. (2026).** *A Permutation-Based Intensity Test for Comparing Spatial Point Processes: Motorcycle Crash Locations and Perceived Road Risks on a Linear Network.*

---

## Overview for Geographers and Applied Researchers

Comparing two spatial point patterns (for example, **objectively recorded traffic crashes** vs. **hazardous locations subjectively reported by drivers/riders**, or **recorded crime** vs. **perceived fear of crime**) is a common task in geography. 

Traditional methods (such as the Spatial Point Pattern Test – SPPT, Fisher's exact test, or grid-based bootstrap tests) divide the study area into arbitrary polygons or grid squares. This makes them sensitive to the **Modifiable Areal Unit Problem (MAUP)** and causes them to lose statistical power on fine grids or fail when the two datasets have very different sample sizes ($n_1 \neq n_2$).

This repository provides two easy-to-use R functions that overcome these limitations:
1. **Continuous Kernel Smoothing:** Instead of counting points in arbitrary grid squares, the method compares continuous normalized **Kernel Density Estimation (KDE)** surfaces.
2. **Exact Random-Relabeling Permutations:** It pools both datasets together and randomly shuffles their group labels to simulate the null hypothesis ($H_0: f_1(u) = f_2(u)$), which guarantees valid Type I error control even when sample sizes are heavily unbalanced ($n_1 \neq n_2$).
3. **Global False Discovery Rate (FDR) Control:** Using the `GET` package, it constructs a two-sided Global FDR envelope across the entire map or road network, pinpointing exactly *where* the two spatial patterns significantly diverge.

---

## Repository Contents

* **`Comparison_of_2pp_2D.R`** – Contains `test_permutation_kde_env()`, the function for comparing two planar (2D) point patterns (`ppp` objects).
* **`Comparison_of_2pp_LN.R`** – Contains `test_permutation_kde_env_linnet()`, the parallelized function for comparing two point patterns constrained to a linear network such as roads or rivers (`lpp` objects).
* **`Example_Comparison_of_2pp.R`** – A self-contained, reproducible tutorial script. It simulates a $10 \times 10\text{ km}$ road network with objective crashes ($n_1 = 600$) and perceived road risks ($n_2 = 300$), runs both the 2D planar test and the linear network test, and plots the results.

---

## Step 1: Installation of Required R Packages

Before running the scripts, open R or RStudio and install the required packages by copying and pasting the following command into the R Console:

```r
install.packages(c("spatstat", "spatstat.linnet", "GET", "foreach", "doParallel", "doRNG", "sf"))
```

*(Note: The `sf` package is only needed if you load your own GIS shapefiles or GeoPackages.)*

---

## Step 2: Running the Simulated Example (`Example.R`)

To see how the methods work without needing any external data files:

1. Download **`Comparison_of_2pp_2D.R`**, **`Comparison_of_2pp_LN.R`**, and **`Example_Comparison_of_2pp.R`** into the same folder on your computer.
2. Open **`Example_Comparison_of_2pp.R`** in RStudio.
3. Set your working directory to that folder (in RStudio top menu: **Session -> Set Working Directory -> To Source File Location**).
4. Run the entire **`Example_Comparison_of_2pp.R`** script.

The script will simulate two point patterns on a $10 \times 10\text{ km}$ road network with two planted anomalies (an *unrecognized crash hotspot* in the North-West and an *overestimated fear hotspot* in the South-East) and display the resulting significance maps.

---

## Step 3: How to Use the Functions on Your Own GIS Data

> **Important Note on Coordinate Systems (CRS):**  
> Always project your spatial data into a **planar metric coordinate system** (in meters, such as UTM or S-JTSK `EPSG:5514` for Czechia) before running the tests. Do **not** use unprojected GPS coordinates (longitude/latitude in degrees, `EPSG:4326`), because the smoothing bandwidth `h` is specified in map units (meters).

### A. Planar (2D) Comparison (`Comparison_of_2pp_2D.R`)

Use this test when your events can occur anywhere within a 2D study region (polygon) or for a macro-scale regional overview:

```r
source("Comparison_of_2pp_2D.R")

# X1 and X2 must be 'ppp' point pattern objects (from the 'spatstat' package)
# defined on the same observation window W (in meters).
res_2d <- test_permutation_kde_env(
  X1 = X1_crashes,          # First point pattern (e.g., actual crashes)
  X2 = X2_perceived,        # Second point pattern (e.g., perceived risks)
  h = 3000,                 # Smoothing bandwidth in meters (if NULL, estimated automatically)
  n_permutations = 1999,    # Number of permutations (use 1999+ for standard grids)
  fdr_level = 0.05,         # 5% False Discovery Rate level
  dimyx = c(64, 64)         # Evaluation grid resolution (64 x 64 pixels)
)

# Plot the results
plot(res_2d, what = c("obs", "lo", "hi", "lo.sign", "hi.sign"))
```

### B. Linear Network Comparison (`Comparison_of_2pp_LN.R`)

Use this test when your events are strictly constrained to a road, street, or river network (`linnet` / `lpp` objects in `spatstat.linnet`):

```r
source("Comparison_of_2pp_LN.R")

# LPP1 and LPP2 must be 'lpp' objects on the same linear network L (in meters).
res_net <- test_permutation_kde_env_linnet(
  LPP1 = LPP_crashes,       # First network point pattern (e.g., actual crashes)
  LPP2 = LPP_perceived,     # Second network point pattern (e.g., perceived risks)
  h = 1250,                 # Smoothing bandwidth in meters (e.g., 1250 m)
  n_permutations = 1999,    # Number of permutations (use ~2 * number of network pixels for fine national grids)
  fdr_level = 0.05,         # 5% False Discovery Rate level
  dimyx = c(256, 256),      # Grid resolution intersecting the network
  num_cores = NULL          # NULL automatically uses all available CPU cores minus 1
)

# Plot the network results
plot(res_net, what = c("obs", "lo", "hi", "lo.sign", "hi.sign"))
```

---

## How to Interpret the Output Plots

Calling `plot(res, what = c("obs", "lo", "hi", "lo.sign", "hi.sign"))` produces panels showing:

* **`obs`**: The observed difference between the two normalized spatial densities, $T_0(u) = \hat{f}_1(u) - \hat{f}_2(u)$.
* **`lo` and `hi`**: The lower and upper bounds of the 95% Global FDR envelope simulated under the null hypothesis of equal spatial distributions.
* **Red areas (`hi.sign`, where `obs > hi`)**: Locations where the relative density of **Pattern 1 ($X_1$) is significantly higher than Pattern 2 ($X_2$)**. In our road-safety study, these are **unrecognized hazards**—places where actual crashes happen significantly more often than riders expect.
* **Green / Blue areas (`lo.sign`, where `obs < lo`)**: Locations where the relative density of **Pattern 2 ($X_2$) is significantly higher than Pattern 1 ($X_1$)**. In our road-safety study, these are places of **heightened subjective concern**—locations where riders perceive high danger, but historical crash counts are comparatively low.

---

## Citation & Contact

If you use this code in your research, please cite our paper and the underlying `GET` package:

* **Mrkvička, T., & Kraft, S. (2026).** A Permutation-Based Intensity Test for Comparing Spatial Point Processes: Motorcycle Crash Locations and Perceived Road Risks on a Linear Network.
* **Mrkvička, T., & Myllymäki, M. (2023).** False discovery rate envelopes. *Statistics and Computing*, 33(5), Article 109. https://doi.org/10.1007/s11222-023-10275-7
* **Myllymäki, M., Mrkvička, T., Grabarnik, P., Seijo, H., & Hahn, U. (2017).** Global envelope tests for spatial processes. *Journal of the Royal Statistical Society: Series B (Statistical Methodology)*, 79(2), 381–404.
