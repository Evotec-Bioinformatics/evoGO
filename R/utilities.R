# Format text output
trimToWidth <- function(x, width) {
  toTrim <- nchar(x) > width
  x[toTrim] <- paste0(strtrim(x[toTrim], width), "...")
  x
}

# Make a copy of evoGO environment to use it for analysis
copyEnvironment <- function(e) {
  res <- new.env(hash = TRUE, parent = emptyenv())
  objects <- as.list(e, all.names = TRUE)
  for (i in names(objects)) {
    assign(i, objects[[i]], res)
  }
  res
}

# Check for valid character parameters of length one
isValidString <- function(s) {
  !is.null(s) && length(s) == 1 && !is.na(s) && is.character(s) && nchar(s) > 0
}

# Check GO release
validateGoRelease <- function(g) {
  assertthat::assert_that(length(g) == 1 && !any(is.na(as.Date(g, format = "%Y-%m-%d"))),
    msg = "Argument 'goRelease' is expected to be a a character vector with length of one in the format \"2023-04-01\""
  )
  assertthat::assert_that(as.Date(g) >= as.Date("2013-04-01"),
    msg = "GO consortium releases prior to 2013-04-01 are not supported by evoGO"
  )
  as.character(as.Date(g, format = "%Y-%m-%d"))
}

# Check Ensembl release
isValidEnsemblRelease <- function(e) {
  assertthat::assert_that(is.null(e) || (length(e) == 1 && (is.numeric(e))),
    msg = "Argument 'ensemblRelease' is expected to be an integer"
  )
}
