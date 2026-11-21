# =============================================================================
# run_all.R - runs every script in Ch-2 ... Ch-9 and saves what each one prints
#
# How to use:
#   1. Unzip Ch-2.zip ... Ch-9.zip into the same folder as this file, so the folder
#      contains Ch-2, Ch-3, ..., Ch-9 side by side.
#   2. In RStudio or VS Code run:  source("path/to/run_all.R")
#      (or setwd() to that folder and run source("run_all.R"))
#   3. It takes a few minutes. When it finishes, attach run_output.zip in the chat.
#
# What it does:
#   - installs the packages the scripts need that you don't have yet
#     (lmtest, quantreg, leaps, pls)
#   - runs each script in its own fresh R process (Rscript), so nothing in your
#     current workspace is touched and no script depends on another
#   - runs one expression at a time, so an error is recorded and the rest of the
#     script still runs
#   - writes the code and output of each script to run_output/<Ch>_<script>.txt
#     and its plots to run_output/<Ch>_<script>.pdf
# =============================================================================

# The small program each script is run through (written to a temporary file below)
# local({ ... }) keeps the runner's own variables out of the workspace the scripts use
.inner_runner <- r"---(
local({
args <- commandArgs(trailingOnly = TRUE)
f <- args[1]; log_file <- args[2]; pdf_file <- args[3]
options(warn = 1, width = 100)
pdf(pdf_file)
con <- file(log_file, open = "wt")
sink(con); sink(con, type = "message")
n_err <- 0
t0 <- Sys.time()
exprs <- tryCatch(parse(f, keep.source = TRUE), error = function(e) {
  cat("PARSE ERROR:", conditionMessage(e), "\n"); NULL })
if (is.null(exprs)) n_err <- 1
srcs <- attr(exprs, "srcref")
for (i in seq_along(exprs)) {
  code <- paste(as.character(srcs[[i]]), collapse = "\n+ ")
  cat("\n> ", code, "\n", sep = "")
  res <- tryCatch(withVisible(eval(exprs[[i]], envir = globalenv())),
                  error = function(e) {
                    n_err <<- n_err + 1
                    cat("ERROR:", conditionMessage(e), "\n")
                    NULL
                  })
  if (!is.null(res) && res$visible) print(res$value)
}
secs <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
cat("\n#### FINISHED:", n_err, "error(s),", secs, "seconds\n")
sink(type = "message"); sink(); close(con)
while (dev.cur() > 1) dev.off()
})
)---"

.run_all <- function() {
  # Folder that contains this file (falls back to the working directory)
  here <- tryCatch({
    p <- sys.frame(1)$ofile
    if (is.null(p)) getwd() else dirname(normalizePath(p))
  }, error = function(e) getwd())
  if (!dir.exists(file.path(here, "Ch-2"))) here <- getwd()
  if (!dir.exists(file.path(here, "Ch-2")))
    stop("Can't find the Ch-2 ... Ch-9 folders. Put run_all.R next to them, or setwd() there.")

  # Install missing packages once
  needed <- c("faraway", "lmtest", "quantreg", "leaps", "pls")
  missing_pk <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pk) > 0) install.packages(missing_pk)

  out_dir <- file.path(here, "run_output")
  dir.create(out_dir, showWarnings = FALSE)
  runner <- tempfile(fileext = ".R")
  writeLines(.inner_runner, runner)
  rscript <- file.path(R.home("bin"), "Rscript")
  summary_lines <- character(0)

  for (ch in sprintf("Ch-%d", 2:9)) {
    files <- list.files(file.path(here, ch), pattern = "\\.R$", full.names = TRUE)
    for (f in files) {
      stem <- paste0(ch, "_", sub("\\.R$", "", basename(f)))
      log_file <- file.path(out_dir, paste0(stem, ".txt"))
      message("running ", stem, " ...")
      system2(rscript, shQuote(c(runner, f, log_file, file.path(out_dir, paste0(stem, ".pdf")))))
      # Read the last line of the log to report errors and timing
      tail_line <- if (file.exists(log_file)) tail(readLines(log_file, warn = FALSE), 1) else "no log"
      line <- sprintf("%-45s %s", stem, tail_line)
      summary_lines <- c(summary_lines, line)
      message("   ", tail_line)
    }
  }

  # Versions of R and the packages, then zip everything for attaching
  writeLines(c(summary_lines, "", capture.output(sessionInfo()),
               paste("faraway version:", as.character(packageVersion("faraway"))),
               paste("pls version:", as.character(packageVersion("pls")))),
             file.path(out_dir, "00_summary.txt"))
  zip_file <- file.path(here, "run_output.zip")
  ok <- tryCatch({
    utils::zip(zip_file, list.files(out_dir, full.names = TRUE), flags = "-j -q")
    file.exists(zip_file)
  }, error = function(e) FALSE, warning = function(w) FALSE)
  if (ok) {
    message("\nDone. Attach this file in the chat: ", zip_file)
  } else {
    message("\nDone. Zip the folder ", out_dir, " and attach it in the chat.")
  }
}

.run_all()
