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
h5adMatrixLoad <- function(path, sparse = TRUE) {
  data <- h5read(path, "X", read.attributes = TRUE)
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

# You can pass the anndata obs or var to this as data to get a data.frame back
# It won't try and do any data conversion, so you might have with chr dtypes
anndataDataframe <- function(data) {
  index <- data$`_index`

  processLists <- function(item) {
    ret <- if(is.list(item) && all(c("categories", "codes") %in% names(item))) {
      colFactor <- factor(item$codes)
      levels(colFactor) <- item$categories
      colFactor
    } else {
      as.numeric(item)
    }
    return(ret)
  }
  # Convert relevant data to factors
  listdata <- lapply(data, processLists)
  dataframe <- data.frame(Reduce(cbind, listdata), row.names = index)
  colnames(dataframe) <- names(listdata)
  # remove index
  dataframe$`_index` <- NULL
  # Remove redundant x and y if present
  dataframe$x <- NULL
  dataframe$y <- NULL
  return(dataframe)
}
