#' @include generics.R
#' @include utilities.R
#' @import flowCore
#' @import magrittr
#' @importFrom rlang .data
#'
#'
NULL


#'
#' UMAP Dimension Reduction
#'
#' Runs the Uniform Manifold Approximation and Projection (UMAP) dimension reduction technique on the flow cytometry data.
#'
#' @param flow_object A Flow Object
#' @param transformation Data transformation method to apply to the data. Possible choices are 'asinh', 'none', 'log' or 'cytofAsinh'. See `details` for more information
#' @param n_neighbors Number of neighbors. See documentation from the `uwot::umap()` function for more details.
#' @param scale Logical argument to determine whether or not to scale the data prior to umap. Defaults to TRUE
#' @param dims Number of dimensions to reduce the data unto. 2 is the best for the vast majority of cases, but 3 can be useful for trajectory analysis.
#'
#' @details
#' Possible choices for data transformation are hyperbolic sine (asinh)
#' using cofactors inferred by `find_cofactors()`, biexponential (none), logarithmic (in base 10), or CYTOF asinh (cofactor of 5).
#'
#' @return A Flow Object
#'
#' @export
#'


UMAP_flow <- function(flow_object,
                      transformation = "asinh",
                      n_neighbors = 15,
                      scale = TRUE,
                      dims = 2) {

  coll <- checkmate::makeAssertCollection()
  if(methods::is(flow_object) != "flow_object") { coll$push("flow_object not recognized. Please supply a valid flow object. ")}
  checkmate::assertCharacter(transformation, len = 1, null.ok = F, any.missing = F, .var.name = "transformation", add = coll)
  checkmate::assertNumeric(n_neighbors, len = 1, lower = 3, upper = 50, null.ok = F, any.missing = F, .var.name = "n_neighbors", add = coll)
  checkmate::assertNumeric(dims, len = 1, lower = 2, upper = 3, null.ok = F, any.missing = F, .var.name = "dims", add = coll)
  checkmate::assertLogical(scale, len = 1, null.ok = F, any.missing = F, .var.name = "scale", add = coll)
  checkmate::reportAssertions(coll)

  out_object <- flow_object
  flowDF <- transform_flow_data(flow_object, transform = transformation)

  if(scale == TRUE){
    flowDF <- scale(flowDF, scale = T, center = T)
  }

  print("Performing dimension reduction using UMAP...")
  umap_mat <- as.data.frame(uwot::umap(flowDF, n_neighbors = n_neighbors, n_components = dims)) %>%
    dplyr::rename(UMAP1 = 1,
                  UMAP2 = 2)
  row.names(umap_mat) <- row.names(flowDF)


  out_object$dims <- umap_mat

  out_object$parameters$transformation <- transformation
  return(out_object)

}

#'
#' Phenograph Clustering
#'
#' Runs Phenograph clustering on the flow data to identify louvain clusters from a flow object with embedded UMAP coordinates.
#' Uses the same transformation used for `UMAP_flow()`.
#'
#' @param flow_object A Flow Object
#' @param k number of neighbors for Phenograph. Defaults to 60.
#' @param scale Logical argument to determine whether or not to scale the data prior to umap. Defaults to TRUE
#'
#' @return A Flow Object
#'
#' @export


Cluster_flow <- function(flow_object,
                         k = 60,
                         scale = TRUE){

  coll <- checkmate::makeAssertCollection()
  if(methods::is(flow_object) != "flow_object") { coll$push("flow_object not recognized. Please supply a valid flow panel. ")}
  checkmate::assertNumeric(k, lower = 5, upper = 150, any.missing = F, len = 1, null.ok = F, .var.name = "k", add = coll)
  checkmate::assertLogical(scale, len = 1, null.ok = F, any.missing = F, .var.name = "scale", add = coll)
  checkmate::reportAssertions(coll)

  out_object <- flow_object
  #flowDF_annot <- get_flowSet_matrix(flow_object, add_sample_id = T)

  if(length(intersect(c("UMAP1", "UMAP2"), colnames(get_flowSet_matrix(flow_object, add_sample_id = T)))) < 2){
    coll$push(paste("UMAP model not found in flow object. First run ",sQuote("UMAP_flow"), ". ", sep = ""))
  }
  checkmate::reportAssertions(coll)

  print("Performing clustering using RPhenograph...")

  transformation <- out_object$parameters$transformation
  flowDF <- transform_flow_data(flow_object, transform = transformation)

  if(scale == TRUE){

    flowDF <- scale(flowDF, scale = T, center = T)
  }


  louvain_df <- data.frame("louvain" = as.character(igraph::membership(Rphenograph::Rphenograph(flowDF, k = k)[[2]]))) %>%
    dplyr::mutate(partition = "1")
  row.names(louvain_df) <- row.names(flowDF)

  cluster_annot <- data.frame("louvain" = as.character(sort(unique(as.numeric(louvain_df$louvain))))) %>%
    dplyr::ungroup()


  out_object$louvain <- louvain_df
  out_object$parameters$Phenograph_k <- k
  out_object$parameters$louvain_annotations <- cluster_annot
  return(out_object)

}


