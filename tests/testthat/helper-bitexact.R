# Some fixtures in this suite are frozen bit for bit: md5 digests of a whole
# fit (helper-irls-driver.R) and hex float literals from the dispersion MLE.
# They are the gate that the compiled core still reproduces the R it replaced,
# and they have caught real regressions -- the Z'WZ triangle, and allowing FMA
# contraction in the linear predictor.
#
# They cannot hold on a machine other than the one that froze them. exp(), log()
# and lgamma() are not bit-identical across libm implementations, so a fit
# drifts in its last bits, and a last-bit difference can cross a convergence
# threshold: on Bioconductor's Linux builders the breast cancer fit stops at a
# different iteration and even the n_iter digest moves. Suppressing fused
# multiply-add (see ./configure) removes the compiler's contribution and makes
# every comparison against an R reference agree, but it cannot make two C
# libraries agree.
#
# So these fixtures run where they were frozen and skip elsewhere, saying so.
# The fingerprint includes the R version because a libm can change under it.
bitexact_fingerprint <- function() {
  paste(R.version$platform, paste0(R.version$major, ".", R.version$minor))
}

# The machine the fixtures in this suite were frozen on. Update it in the same
# commit that regenerates them, never on its own.
bitexact_reference <- "aarch64-apple-darwin20 4.5.0"

skip_unless_bitexact_platform <- function() {
  testthat::skip_if(
    !identical(bitexact_fingerprint(), bitexact_reference),
    sprintf("bit-exact fixtures were frozen on %s; this is %s",
            bitexact_reference, bitexact_fingerprint()))
}
