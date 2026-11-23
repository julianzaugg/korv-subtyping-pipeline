# Tab-separated file reading and writing on top of data.table.
#
# readr/vroom is deliberately avoided: its compiled library crashes on loading on some platforms
# (Apple Silicon conda builds), and data.table is faster on the multi-million line vsearch files.
# The functions below keep the small subset of the readr interface the pipeline scripts use.
suppressPackageStartupMessages({
  library(data.table)
  library(tibble)
})

cols <- function(...) list(...)
cols_only <- function(...) structure(list(...), only = TRUE)

.class_of <- function(code) c(c = "character", d = "numeric", i = "integer")[[code]]

# col_names: TRUE (header in file) or a character vector of names (no header in file)
# col_types: NULL, a string such as "ccd", cols(.default = "c") or cols_only(X1 = "c", X9 = "c")
read_tsv <- function(file, col_names = TRUE, col_types = NULL, progress = FALSE, ...) {
  header <- isTRUE(col_names)
  names_given <- if (is.character(col_names)) col_names else NULL
  spec_only <- is.list(col_types) && isTRUE(attr(col_types, "only"))

  if (file.size(file) == 0) {   # an empty file, e.g. no search hits
    n <- if (!is.null(names_given)) names_given else character(0)
    cls <- if (is.character(col_types) && length(col_types) == 1) vapply(strsplit(col_types, "")[[1]], .class_of, "") else rep("character", length(n))
    return(as_tibble(setNames(lapply(cls, function(k) vector(k, 0)), n)))
  }

  args <- list(file = file, sep = "\t", header = header, quote = "", na.strings = c("", "NA"),
               data.table = FALSE, showProgress = FALSE, check.names = FALSE)
  if (is.character(col_types) && length(col_types) == 1) {
    args$colClasses <- unname(vapply(strsplit(col_types, "")[[1]], .class_of, ""))
  } else if (is.list(col_types) && !spec_only) {
    args$colClasses <- "character"
  }
  if (spec_only) {
    keep <- names(col_types)
    idx <- as.integer(sub("^X", "", keep))        # unnamed columns are X1, X2, ...
    args$select <- idx
    args$colClasses <- "character"
    out <- do.call(fread, args)
    names(out) <- keep
    return(as_tibble(out))
  }
  out <- do.call(fread, args)
  if (!is.null(names_given)) names(out) <- names_given
  as_tibble(out)
}

write_tsv <- function(x, file, na = "NA", ...) {
  fwrite(as.data.frame(x), file, sep = "\t", quote = FALSE, na = na, compress = "auto")
  invisible(x)
}
