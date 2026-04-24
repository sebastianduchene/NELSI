suppressMessages({
  library(testthat)
  library(ape)
  library(castor)
  library(phangorn)
  library(geiger)
})

# Source all NELSI R files
rfiles <- list.files("R", full.names = TRUE, pattern = "\\.R$")
for (f in rfiles) source(f)

# Run tests
test_results <- testthat::test_file(
  "tests/testthat/test-castor-integration.R",
  reporter = "progress"
)
