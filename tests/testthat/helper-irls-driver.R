# A bit-exact digest of a fit's loop outputs.
#
# The gate needs to detect ANY change in these matrices, and it used to do that
# by fitting the same data through the old R loop and comparing. That loop is
# gone -- it was reachable only under fuse_rho, which existed only to paper over
# the loop being in R -- so the oracle is now a frozen digest instead.
#
# md5 of serialize(), not sum() or a tolerance: U, se_U and tau_g_array are
# q x genes, so storing them whole cost 960 KB of package for one fixture, while
# a summary statistic would let a real change hide behind a collision. serialize()
# with the default xdr = TRUE writes big-endian IEEE doubles, so the digest is
# the same on every platform R builds on.
irls_digest <- function(fit) {
  fields <- c("B", "U", "se_B", "se_U", "alpha", "tau_g_array", "tau_blocks",
              "percell_bleed_rho", "history", "n_iter", "converged")
  vapply(fields, function(field) {
    path <- tempfile()
    on.exit(unlink(path), add = TRUE)
    connection <- file(path, open = "wb")
    serialize(fit[[field]], connection, xdr = TRUE)
    close(connection)
    unname(tools::md5sum(path))
  }, character(1))
}
