# PACE 0.99.3

* The portability tolerance from 0.99.2 was set at 1e-12, which is below the
  noise floor of the quantity it gates: these are sums over ~7,900 cells, and
  n * eps is 1.7e-12. aarch64 failed at 2.8e-12. Raised to 1e-10, the gate the
  rest of the suite already uses.

# PACE 0.99.2

Portability fixes found by Bioconductor's multi-platform builders. No change to
any reported quantity: breast cancer holds at 1,637 calls, melanoma at 46 with
SPP1 Macrophage<-Tumour at -0.066840.

* Fused multiply-add is now suppressed by a compiler flag rather than by
  pragmas. `#pragma GCC optimize("fp-contract=off")` is a documented debugging
  aid that does not reliably override GCC's default, and `gene_solve.cpp` --
  the solver -- never included the header at all, so contraction was live there
  even on macOS. `./configure` probes the compiler and writes `-ffp-contract=off`
  into `src/Makevars`, leaving it empty where the flag is rejected. This took
  all six platforms from failing to x86_64 passing.

* Comparisons between the compiled core and the R it replaces are made at a
  tolerance rather than bit for bit where the quantity is a REDUCTION -- the
  kernels, the ambient field, the linear predictor's Z sum. A reduction's
  accumulation order is the compiler's to choose, and AVX and NEON choose
  differently; aarch64 failed at 1e-15 where x86_64 passed. Shapes, names and
  which entries are zero must still match exactly.

* The frozen fixtures -- the IRLS digests and the pinned dispersion values --
  run on the machine that froze them and skip elsewhere. They compare against
  stored constants, and `exp`, `log` and `lgamma` are not bit-identical across
  libm implementations, so no compiler setting can make them portable.

* The shrinkage tests pin `shrink_threads = 1`. The default is 4, and on a
  platform without `fork()` the shrinkage says so in a message that arrived
  before the one the tests assert on, failing every Windows build.

# PACE 0.99.1

Correctness fixes from an internal audit of the fitting engine. All of these
produced wrong numbers with no error raised. The analyses in the accompanying
manuscript were produced with the pre-fix engine and were checked to be
unaffected.

Reported output is unchanged. On the shipped breast cancer subset the shrunken
slopes move by at most `1.6e-06` and the local false sign rates by `3.2e-05`
across all 17,792 (gene, focal, neighbour) rows, with the same 60 calls at
`lfsr < 0.05`, no flipped calls, and no sign changes. The raw solver output
moves at the `1e-04` level, and moves *towards* the full-double reference
(fixed effects `8.4e-04 -> 5.3e-05` from it), because the ridge no longer
drops into single precision. Every other fix is exactly reproducing: with the
precision guard disabled the refit is identical to 0.99.0 (`max|diff| = 0`).

## Correctness fixes

* An `NA` in any random-effect variable no longer corrupts the design matrix.
  `model.matrix()` drops such rows, which left the sparse assembly recycling
  covariates so that cells received other cells' values; this now stops with a
  message naming the column and the number of rows dropped.
* The condition interaction term is derived from `condition_col` instead of a
  hardcoded `"Responder"` prefix. Cohorts whose condition column had any other
  name silently produced an all-zero condition indicator, and therefore an
  identically zero responder spatial block, while the output looked complete.
  A `resp_term` that names no level of `condition_col` is now refused.
* Standard errors are computed on every iteration, so a fit that exits early at
  convergence no longer returns all-`NA` `se_B`/`se_U` (which the shrinkage
  layer then discarded). Convergence on a float interior takes one further pass
  in double so the retained standard errors are not float-quantised.
* Non-finite genes are counted and reported instead of being masked out of the
  convergence metric, and a fit that has lost genes is no longer reported as
  converged. In the streaming engine the NaN repair is re-checked after it runs,
  and genes that remain non-finite are reported on every iteration, including
  the last, where the interior is already double and the repair cannot help.
* The ridge is kept out of single precision. It is accumulated inside the float
  kernel, where the smallest representable increment on a diagonal of order
  `sum(w)` is `eps_float * sum(w)`; below that the regularisation is quantised
  away, at a threshold that falls as `1/n`, so the same model and design lose
  their ridge purely by growing the cohort. The interior is now forced to double
  whenever the ridge is small relative to the diagonal scale.
* The variance components are bounded from above by a new `tau_max` argument
  (default `100`). Every previous tau write was a `pmax` against a floor with no
  ceiling, so a column identified only by the ridge climbed by `mean(u^2)` per
  iteration while its inflated posterior variance fed the next update. The
  default sits roughly 6-30x above the largest component observed on healthy
  cohorts and just below `var(z)`, the variance of the working response that the
  components decompose; it does not bind on any fit in the test suite. A binding
  cap is reported, since it means the term is not identified by the data.
* The empirical-Bayes floor in the adaptive tau update is scaled by the median
  across genes rather than the mean, so a single runaway gene can no longer lift
  the floor for every gene in its row. The shrinkage target itself is unchanged.
* The per-iteration largest variance component is recorded in the fit history and
  printed alongside the block medians, since a runaway confined to one term does
  not move a median taken over the whole block.
* mash is now given the correlation between conditions under the null. The
  focal cell types are conditions estimated in the same per-gene solve, so
  their errors are correlated, and the previous default of `V = I` treated them
  as independent and mis-calibrated the local false sign rates. The estimate is
  taken over the conditions that vary: sparse focal-neighbour pairs are zeroed
  by the kernel drop without their column being removed, so those conditions
  have a constant z column, and `cov2cor()` on the full matrix returns `NaN`
  for every entry. Degenerate conditions are left uncorrelated, which is the
  correct limit, and the estimate falls back to independence when fewer than
  two conditions vary. Set `null_correlation = FALSE` in `paceShrink()` to
  restore the old behaviour.

  On the shipped breast cancer subset the estimated correlations reach 0.20 in
  absolute value, the shrunken slopes correlate at 0.9954 with the independent
  fit, and one call changes net at `lfsr < 0.05` (60 to 61: CD163 in macrophages
  near tumour and PECAM1 in B cells near T cells gained, ERBB2 in tumour near
  T cells lost, all three within 0.02 of the threshold). The 1,294 nominal sign
  changes are numerical noise: none is significant in either pass and their
  median absolute effect is 4e-21.
* Leave-one-patient-out variance decomposition now subsets the contamination
  offset matrix alongside `mu`. Leaving it at full length offset every cell
  after the dropped patient's first row, so each leave-one-out row was wrong
  while the full-data row was correct.
* The condition-by-proximity pair contributions are matched to the neighbours
  whose interaction terms are actually present, rather than to the first
  `n` entries of the neighbour list, and pair rows are labelled from the same
  matching vector. The previous indexing misaligned neighbour columns whenever
  an interaction term was missing.
* `plotResponseCurve()` draws each arm's fitted slope on the arm it belongs to.
  It took the arm carrying the responder interaction from the first level of
  the plotted factor, which is the alphabetically first one, while the model's
  reference is the level recorded in `params$resp_term`. The two agree only by
  coincidence, so on a cohort whose arms are `"R"` and `"NR"`, where the model
  references `R` but `NR` sorts first, both dashed lines were drawn on, and
  labelled with, the wrong arm. The agreement was not even stable for one
  cohort: `"PD"` sorts before `"nonPD"` under `C` collation and after it under
  `en_US`, so the same call could label the figure differently on two machines,
  and the `condition` vignette exercises exactly that pair. A `resp_term` that
  does not name a level of `condition_col` is now an error rather than a silent
  guess. No fitted quantity is affected: the binned means are computed from the
  data, and the slopes, decomposition and pair contributions all index the
  coefficient matrices by explicit term name.

## New features

* The memory and approximation settings the fitter accepts are documented on
  `paceModel()`, under "Memory and approximation settings". Three of them were
  reachable but undocumented: `ambient_mode`, which decides whether the ambient
  field is materialised or recomputed per chunk and is what lets a full
  transcriptome panel fit in memory; `alpha_warmup`, which freezes the
  dispersion MLE after its first few iterations and **defaults to 6**, so every
  fit already runs that approximation; and `alpha_max_n`, which caps the cells
  the dispersion MLE sees. The note on `alpha_max_n` records that a cap
  validated on two cohorts of roughly 10^5 cells did not hold on a 1.2M cell
  panel, where it left the number of calls unchanged while changing which
  calls they were.

* `cellContamination()` reports, per cell, the contamination loading `rho_i` and
  the fraction of the cell's expected counts attributed to its local ambient
  field, `sum_g mu_spill_ig / sum_g mu_ig`. This is the quantity the solver
  already printed in its fitting trace and then discarded; both engines now
  retain it as `contam_frac`, which costs two length-`n` accumulators rather
  than the `n x G` matrices, so it is available from an ordinary streaming fit
  without `return_mu`. A high fraction marks a cell whose profile is largely
  explained by its neighbours, the expected signature of a segmentation or
  transcript-assignment error.

  Read the fraction rather than the loading. The ambient field is
  cross-cell-type, so a cell with no differently-typed neighbour inside the
  technical kernel has no ambient signal at all; its `rho_i` is unidentified and
  shrinks to the empirical Bayes prior mean, giving every such cell the same
  apparently middling loading. On the shipped breast cancer subset that is 4,589
  of 7,898 cells, including 84% of the tumour cells.

* `paceDecompose()` rebuilds the fitted means it needs instead of requiring them
  to have been stored. `mu` is a deterministic function of what the fit already
  holds: on the log scale `eta = X B + Z U` has rank at most `p + q`, so `B`,
  `U`, `X_fixed` and `re_meta$Z` are its factored form and the dense `n x G`
  matrix is the expanded copy. Only the ambient field is recomputed, from the
  counts in the `SpatialExperiment` the function already takes. The rebuild is
  exact rather than approximate: `max|diff| = 0` against the solver's own
  matrices, and a decomposition identical to the stored one, across per-cell
  contamination and no contamination, edge correction on and off, a non-default
  technical bandwidth, per-image kernels, a condition cohort, and a non-
  alphabetical cell-type order.

  `technical_offset_mat` is rebuilt alongside `mu`, not as a detail: it gates
  the spillover block, so rebuilding only `mu` leaves a decomposition that
  reports no spillover without raising anything. On the packaged fit that is a
  block with a median of 1.6% and a maximum of 46%.

  The ambient field is rebuilt from the settings the fit recorded, never from
  the defaults of the exported `ambientField()`, since a fit made with another
  bandwidth, image grouping, cell-type order or edge correction would otherwise
  be handed a different field and return plausible, wrong numbers. `paceModel()`
  now records `edge_correct`, the one such input that was not already kept; a
  fit made before this is refused rather than rebuilt on a guess.

  A fit can therefore be saved without its `n x G` matrices and still be
  re-decomposed. Dropping them takes the packaged example fit from 59.8 MB to
  7.1 MB in memory and 1.03 MB on disk. Note that this saves storage, not peak
  memory: the matrices are still materialised while the decomposition runs.

* A second vignette covers condition-stratified cohorts: what `condition_col`,
  `kernel_per_image` and `image_re` do, the responder spatial state block in the
  decomposition, `pairVariance()` / `plotPairHeatmap()` with
  `block = "responder"`, and the gene-level responder slopes.

  It runs on a new shipped subset of the CosMx melanoma cohort of Dong et al.
  (doi:10.5281/zenodo.14708000, CC-BY-4.0). The subset keeps every cell and
  every patient of the published cohort, 56,274 cells across 26 sections, and
  reduces only the panel, to the 180 most widely detected genes. File size is
  cells times genes, and of the two it is the cells that have to be preserved:
  the model reads local neighbourhoods, and the condition effect is a
  between-patient contrast. Cropping cells was tried first and abandoned, since
  at 15% of each section nothing reached `lfsr < 0.05` at all and no cell-based
  crop under about 5 MB recovered the result.

  The vignette therefore computes the melanoma finding rather than citing it:
  macrophage *SPP1* responds less steeply to tumour proximity in progressive
  disease, `-0.0663` with an `lfsr` of 0, the second strongest responder call in
  the fit. `inst/scripts/make-mel-cosmx-subset.R` records the derivation.

  The reduced panel leaves estimates alone but shifts mash's calibration, which
  the vignette says: 76 responder calls here against 46 on the full panel, and
  SPP1 at `-0.0663` against `-0.0668`. Below about 150 genes a focal type runs
  out of genes to decompose and `paceDecompose()` fails, so the panel cannot be
  cut much further.

* Three plotting additions for condition cohorts, ported from the code behind
  the manuscript's melanoma figures.

  `plotProximity()` gains a `condition` argument, dodging one box per arm within
  each proximity bin instead of pooling the arms. `TRUE` uses the fit's own
  `condition_col`.

  `plotResponseCurve()` draws the population view: binned means and standard
  errors per arm, with the PACE slopes overlaid as dashed lines read from the
  fit rather than smoothed from the points. Its x axis is the model's own
  covariate, the Gaussian kernel density, so the slope is drawn on the scale it
  was estimated on.

  `plotResponseMap()` gives the spatial view for one exemplar section per arm:
  a neighbour-density kernel map with the focal cells coloured by expression,
  above expression against density with a fitted line. Colour and axis scales
  are shared across arms so the panels are comparable.

# PACE 0.99.0

* Initial package: `paceFit()` for SpatialExperiment input, the `PACEFit` S4
  class, and the `neighbourSlopes()`, `varianceDecomposition()`, and
  `topDrivers()` accessors, wrapping the locked PACE statistical engine.
