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
  !is.null(s) && length(s) == 1 && !is.na(s) && is.character(s)
}
