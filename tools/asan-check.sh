#!/bin/bash
# AddressSanitizer check for the compiled core.
#
# Why this exists: the fixture gates compare numbers from runs that finished, so
# they cannot see a memory error that does not happen to crash. One did exactly
# that -- rows written into an unsized Eigen matrix passed all 20 configurations
# at 2.59e-13 and segfaulted only on a full cohort. -DNDEBUG (from ~/.R/Makevars)
# disables Eigen's own bounds assertions, so nothing else catches it either.
#
# ASan instruments every access, so it finds the fault on the FIRST execution of
# the bad line. The workload is therefore deliberately tiny: a few hundred cells
# and two iterations exercise the solve, the linear predictor, the working
# response and the dispersion path, which is all that is needed. It is also all
# that is affordable -- ASan makes Eigen's inner loops hundreds of times slower,
# and a full fit does not finish in a useful time.
#
# Usage:  tools/asan-check.sh [library-directory]
# Exit:   0 clean, 1 build failure, 2 AddressSanitizer reported an error.
set -u
package_root="$(cd "$(dirname "$0")/.." && pwd)"
library="${1:-${TMPDIR:-/tmp}/pace-asan-lib}"
workload="${TMPDIR:-/tmp}/pace-asan-workload.R"
makevars="${TMPDIR:-/tmp}/pace-asan-makevars"
log="${TMPDIR:-/tmp}/pace-asan.log"

# R's own compiler, so the sanitizer runtime matches the one it would link.
compiler="$(cd "$package_root" && R CMD config CXX17 | awk '{print $1}')"
runtime="$("$compiler" -print-file-name=libclang_rt.asan_osx_dynamic.dylib 2>/dev/null)"
if [ ! -f "$runtime" ]; then
  echo "asan-check: no sanitizer runtime for $compiler; on Linux the runtime is"
  echo "            linked into the object and LD_PRELOAD is not needed -- drop"
  echo "            the DYLD_INSERT_LIBRARIES below and run the same workload."
  exit 1
fi

# Extend the user's Makevars rather than replacing it: R_MAKEVARS_USER overrides
# wholesale, which would silently drop their compiler and OpenMP settings.
{ [ -f "$HOME/.R/Makevars" ] && cat "$HOME/.R/Makevars"
  cat <<'FLAGS'

# ---- appended by tools/asan-check.sh ----
# -O2, not the -O1 ASan documents: -O1 leaves Eigen unvectorised and the fit
# becomes unusably slow. LDFLAGS carries the sanitizer into the link, which
# CXXFLAGS alone does not do for an R shared object.
CXXFLAGS   += -g -O2 -fsanitize=address -fno-omit-frame-pointer
CXX17FLAGS += -g -O2 -fsanitize=address -fno-omit-frame-pointer
LDFLAGS    += -fsanitize=address
FLAGS
} > "$makevars"

mkdir -p "$library"
# Stale objects are not rebuilt just because the flags changed, so clean first.
rm -f "$package_root"/src/*.o "$package_root"/src/core/*.o "$package_root"/src/PACE.so
# --no-test-load: R CMD INSTALL runs its load test through /bin/sh, and macOS
# strips DYLD_* when a SIP-protected binary starts, so ASan is not preloaded and
# the test fails with "interceptors are not working".
if ! R_MAKEVARS_USER="$makevars" R CMD INSTALL --no-multiarch --no-test-load \
     --library="$library" "$package_root" > "$log" 2>&1; then
  echo "asan-check: instrumented build FAILED"; tail -20 "$log"; exit 1
fi

cat > "$workload" <<'WORKLOAD'
suppressPackageStartupMessages(library(PACE, lib.loc = Sys.getenv("PACE_ASAN_LIB")))
spe <- readRDS(system.file("extdata", "bc_xenium_subset.rds", package = "PACE"))
x <- SpatialExperiment::spatialCoords(spe)[, 1]
spe <- spe[, x <= stats::quantile(x, 0.08)]          # a few hundred cells
spe <- spe[seq_len(min(24L, nrow(spe))), ]           # and a handful of genes
set.seed(1)
fit <- paceModel(spe, celltype_col = "cellType", image_col = "sample_id",
                 contamination = "percell_hc", dispersion = "nb1",
                 threads = 2L, chunk_size = 8L, n_iter = 2L, min_iter = 1L,
                 verbose = FALSE)
cat("asan workload: cells", ncol(spe), "genes", nrow(spe),
    "iterations", fit@fit$n_iter, "\n")
WORKLOAD

# Rscript is a real Mach-O binary, so DYLD_INSERT_LIBRARIES survives here even
# though it does not survive R CMD's shell wrapper.
rscript="$(R RHOME)/bin/Rscript"
DYLD_INSERT_LIBRARIES="$runtime" \
ASAN_OPTIONS=detect_leaks=0:detect_container_overflow=0:abort_on_error=0 \
PACE_ASAN_LIB="$library" \
  "$rscript" "$workload" > "$log" 2>&1
status=$?

if grep -q "ERROR: AddressSanitizer" "$log"; then
  echo "asan-check: FAILED -- AddressSanitizer reported an error"
  grep -A25 "ERROR: AddressSanitizer" "$log" | head -30
  exit 2
fi
if [ $status -ne 0 ]; then
  echo "asan-check: workload exited $status without an AddressSanitizer report"
  tail -20 "$log"; exit 2
fi
grep "asan workload:" "$log"
echo "asan-check: clean"
exit 0
