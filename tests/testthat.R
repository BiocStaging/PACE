library(testthat)
library(PACE)

# Windows has failed this suite five builds running with a SILENT process death:
# no error, no testthat tally, no "Execution halted". R installs its SIGSEGV
# handler on Unix-alikes only, so an access violation on Windows kills the
# process with nothing printed, and R CMD check then keeps only the last dozen
# lines of whatever had been flushed.
#
# That truncation actively misleads. Every progress line PACE prints goes
# through cat(), i.e. stdout, which R CMD check block-buffers -- so up to 4 KB is
# lost on an abnormal exit, and the surviving lines named a test that had
# finished long before the crash. stderr is unbuffered, so a line written there
# before each file and each test survives the kill and names what was actually
# running. The built-in progress reporters all report a file on FINISH, which is
# exactly the wrong end for this.
#
# CheckReporter stays for its tally and failure summary. It is not what fails the
# check -- test_check()'s stop_on_failure does that, and R CMD check renames
# testthat.Rout to .fail itself on a non-zero exit. The breadcrumbs go alongside
# it, not instead of it.
BreadcrumbReporter <- R6::R6Class(
  "BreadcrumbReporter",
  inherit = testthat::Reporter,
  public = list(
    start_file = function(filename) {
      cat(sprintf("[breadcrumb] FILE %s\n", filename), file = stderr())
    },
    start_test = function(context, test) {
      cat(sprintf("[breadcrumb]   test: %s\n", test), file = stderr())
    }
  )
)

test_check("PACE",
           reporter = testthat::MultiReporter$new(
             reporters = list(testthat::CheckReporter$new(), BreadcrumbReporter$new())))
