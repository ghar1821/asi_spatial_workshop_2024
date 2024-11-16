library(rhdf5)

# Scipy's CSR format can be read into a dgCMatrix with some slight adjustments
convertCSR <- function(matrix_data) {
  # Extract the components
  data <- matrix_data$data        # Non-zero values
  indices <- matrix_data$indices  # Column indices of non-zero values
  indptr <- matrix_data$indptr    # Row pointers for non-zero values
  shape <- attr(matrix_data, "shape") # shape stored as attribute

  # Reconstruct the CSR matrix (have to add 1 to indices for R)
  csr_matrix <- sparseMatrix(
    i = as.integer(rep(seq_len(length(indptr)-1), diff(indptr))),  # expand row indices
    j = as.integer(indices + 1),  # R uses 1-based indexing, so add 1 to column indices
    x = as.numeric(data),
    dims = shape,
  )
}

# h5ad can have dense matrices already, we'll test if we actually need to convert here
h5adMatrixLoad <- function(path, group = "X", sparse = TRUE) {
  data <- h5read(path, group, read.attributes = TRUE)
  if (is.matrix(data)) {
    # Anndata matrix uses columns for features, opposite of R
    message("Matrix already in dense format")
    if (sparse) {
      message("Converting to sparse matrix (use sparse = FALSE to disable this)")
      data <- Matrix(data, sparse = TRUE)
    }
  } else if (
      is.list(data) && all(c("data", "indices", "indptr") %in% names(data))
  ) {
    message("Matrix is in CSR format, performing conversion")
    data <- t(convertCSR(data))
  } else {
    stop("Unknown matrix format, please handle manually.")
  }
  return(data)
}

anndataDataframe <- function(path, group) {

  data <- h5read(path, group, read.attributes = TRUE)

  if (!"_index" %in% names(data)) {
    warning("No _index found in data, using sequential row names")
    index <- seq_len(length(data[[1]]))
  } else {
    index <- data$`_index`
  }

  processLists <- function(item) {
    if (is.list(item)) {
      if (all(c("categories", "codes") %in% names(item))) {
        # Handle categorical data
        colFactor <- factor(item$codes + 1)
        levels(colFactor) <- item$categories
        return(colFactor)
      } else if ("dtype" %in% names(item)) {
        # Handle typed data
        switch(item$dtype,
               "int32" = as.integer(item$values),
               "int64" = as.integer(item$values),
               "float32" = as.numeric(item$values),
               "float64" = as.numeric(item$values),
               "bool" = as.logical(item$values),
               item$values  # default case
        )
      } else {
        warning("Unknown list structure encountered")
        return(unlist(item))
      }
    } else {
      # Preserve numeric types
      if (is.numeric(item)) {
        storage.mode(item) <- if (all(item == floor(item))) "integer" else "double"
      }
      return(item)
    }
  }

  # Process each column while maintaining structure
  listdata <- lapply(data, processLists)

  # Create empty dataframe with correct number of rows
  n_rows <- length(index)
  dataframe <- data.frame(row.names = index,
                          check.names = FALSE,  # Preserve column names as-is
                          stringsAsFactors = FALSE)  # Don't auto-convert strings to factors

  # Add columns one by one to preserve types and attributes
  excluded_cols <- c("_index", "x", "y")
  for (col_name in names(listdata)) {
    if (!col_name %in% excluded_cols) {
      dataframe[[col_name]] <- listdata[[col_name]]
    }
  }

  # Add attributes for metadata if present
  if (!is.null(attr(data, "h5annotations"))) {
    attr(dataframe, "h5annotations") <- attr(data, "h5annotations")
  }

  return(dataframe)
}
