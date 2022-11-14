#' @include generics.R
#' @import flowCore
#' @import magrittr
#' @importFrom rlang .data
#'
NULL

#' Export a Flow Object to a CDS (cell_data_set).
#'
#' Generates a CDS (Monocle3) Object from a Flow Object.Internal function.
#'
#' @param flow_object A Flow Object
#' @param scale Logical argument to determine whether or not to use scaled data. Default is TRUE .
#'
#' @return A spreadsheet file
#'
#'

flow_object_to_cds <- function(flow_object,
                               scale = TRUE){
  coll <- checkmate::makeAssertCollection()

  if(methods::is(flow_object) != "flow_object"){
    coll$push("flow_object not recognized. Please supply a valid flow object.")
  }
  checkmate::reportAssertions(coll)
  if(is.null(flow_object$dims)){
    coll$push(paste("UMAP model not found in flow object. First run ",sQuote("UMAP_flow"), ". ", sep = ""))
  }
  if(is.null(flow_object$louvain)){
    coll$push(paste("Louvain clusters not found in flow object. First run ",sQuote("Cluster_flow"), ". ", sep = ""))
  }
  checkmate::reportAssertions(coll)


  #markerLS <- flow_object$parameters$include_markers
  flowDF <- get_flowSet_matrix(flow_object = flow_object,
                               add_sample_id = T)

  counts <- transform_flow_data(flow_object, transform = flow_object$parameters$transformation)
  #row.names(counts) <- flowDF$cellID
  if(scale == TRUE){
    counts <- scale(counts, scale = T, center = T)
  }
  umap_mat <- flow_object$dims
  # umap_mat <- as.matrix(umap_mat[row.names(counts),])
  #return(umap_mat)
  #return(counts)
  cds <- suppressWarnings(monocle3::new_cell_data_set(expression_data = t(counts)))
  SingleCellExperiment::reducedDims(cds) <- list(UMAP=as.matrix(umap_mat))
  return(cds)
}


#' Set pseudotime partitions
#'
#' Monocle 3 attempts to detect separate groups of cells that likely fall within the same trajectories. By default, all cells are set up to be part of the same
#' partition (named "1"), but this behavior does not always make biological sense. For instance, when analyzing whole blood PBMCs, a trajectory between monocytes and T or B cells is not expected.
#' This will detect them as 2 separate partitions and will make them to be analyzable in parallel.
#' Partitions can also be manually set by providing a data.frame of partitions for each cell, with a 'cellID' column.
#'
#' @param flow_object A Flow Object
#' @param partitions A partitions data.frame to be provided if automatic is set to FALSE, Must contain a cellID column that match the flow object's cell IDs, and a partition column.
#' @param automatic Logical argument to determine whether or not to to automatically detect partitions using monocle3. Defaults to TRUE (recommended)
#' @param scale Logical argument to determine whether or not to scale the data prior to automatic detection of partitions.
#'
#' @return A Flow Object with partition annotations.
#'
#' @export
#'

SetPartitions <- function(flow_object,
                          partitions,
                          automatic = TRUE,
                          scale = FALSE){
  coll <- checkmate::makeAssertCollection()

  if(methods::is(flow_object) != "flow_object"){
    coll$push("flow_object not recognized. Please supply a valid flow object.")
  }
  checkmate::assertLogical(scale, any.missing = F, null.ok = F, len = 1, .var.name = "scale", add = coll)

  if(is.null(flow_object$louvain)){
    coll$push(paste("Louvain clusters not found in flow object. First run ",sQuote("Cluster_flow"), ". ", sep = ""))
  }
  checkmate::reportAssertions(coll)

  out_object <- flow_object
  if(automatic == TRUE){
    cds <- flow_object_to_cds(flow_object, scale = scale)
    cds <- monocle3::cluster_cells(cds)

    #flow_ob
    part_df <- as.data.frame(cds@clusters$UMAP$partitions) %>%
      tibble::rownames_to_column("cellID") %>%
      dplyr::rename(partition = 2)
    out_object$louvain <- flow_object$louvain %>%
      tibble::rownames_to_column("cellID") %>%
      dplyr::select(.data$cellID, .data$louvain) %>%
      dplyr::left_join(.data, part_df, by = "cellID") %>%
      tibble::column_to_rownames("cellID")


  } else{
    part_df <- as.data.frame(partitions) %>%
      tibble::rownames_to_column("cellID") %>%
      dplyr::rename(partition = 2)

    out_object$louvain <- flow_object$louvain %>%
      tibble::rownames_to_column("cellID") %>%
      dplyr::select(.data$cellID, .data$louvain) %>%
      dplyr::left_join(.data, part_df, by = "cellID") %>%
      tibble::column_to_rownames("cellID")

  }
  return(out_object)

}


#' Export Flow Object to a monocle 3 CDS object.
#'
#' Monocle 3 attempts to detect separate groups of cells that likely fall within the same trajectories. By default, all cells are set up to be part of the same
#' partition (named "1"), but this behavior does not always make biological sense. For instance, when analyzing whole blood PBMCs, a trajectory between monocytes and T or B cells is not expected.
#' This will detect them as 2 separate partitions and will make them to be analyzable in parallel.
#' Partitions can also be manually set by providing a data.frame of partitions for each cell, with a 'cellID' column.
#'
#' @param flow_object A Flow Object
#' @param partition Name of a partition in the Flow Object to perform the trajectory analysis with. Partitions are first set to '1' or using the 'SetPartitions' function.
#' @param downsample Number of events to keep PER CLUSTER. Cluster based downsampling helps to reducing the computational weight of high frequency clusters, especially if using high subsampling numbers at object creation.
#' @param min Minimum number of events to take per cluster.
#' @param scale Logical argument to determine whether or not to scale the data prior to automatic detection of partitions.
#'
#' @return A Flow Object with partition annotations.
#'
#' @export
#'

flow_to_cds <- function(flow_object,
                        partition,
                        downsample = 5000,
                        min = 200,
                        scale = FALSE){

  coll <- checkmate::makeAssertCollection()

  if(methods::is(flow_object) != "flow_object"){
    coll$push("flow_object not recognized. Please supply a valid flow object.")
  }
  checkmate::assertCharacter(partition, any.missing = F, null.ok = F, len = 1, .var.name = "partition", add = coll)
  checkmate::assertLogical(scale, any.missing = F, null.ok = F, len = 1, .var.name = "scale", add = coll)
  checkmate::assertNumeric(downsample, any.missing = F, null.ok = F, len = 1,lower = 500, upper = 100000, .var.name = "downsample", add = coll)
  checkmate::assertNumeric(min, any.missing = F, null.ok = F, len = 1,lower = 100, upper = 100000, .var.name = "min", add = coll)

  checkmate::reportAssertions(coll)
  if(is.null(flow_object$louvain$partition)){
    coll$push("No partitions detected in the flow_object. Run `SetPartitions()` first. ")
  }
  checkmate::reportAssertions(coll)
  if(!is.null(partition)){
    invalid.parts <- setdiff(as.character(partition), as.character(unique(flow_object$louvain$partition)))
    if(length(invalid.parts) > 0){
      coll$push(paste("Invalid partition: ", paste(sQuote(invalid.parts), collapse = "; "),". Valid values are ",
                      paste(sQuote(unique(flow_object$louvain$partition)), collapse= "; ")))
    }
    checkmate::reportAssertions(coll)
    parts <- partition
  } else {
    parts <- unique(flow_object$louvain$partition)
  }
  #return(parts)
  out_object <- flow_object

  if(is.null(out_object$trajectories)){
    out_object$trajectories <- list()
    out_object$louvain$pseudotime <- rep(NA, dim(get_flowSet_matrix(flow_object, add_sample_id = T))[1])
  }

  if(length(unique(flow_object$louvain$partition)) == 1){
    message("Warning: Generating graph on whole dataset.. Consider using SetPartitions if the data is too disjointed.")
    cds <- flow_object_to_cds(flow_object, scale = scale)

  } else{
    cds <- flow_object_to_cds(SubsetFlowObject(flow_object = flow_object, subset = partition == parts), scale = scale)
  }

  cds <- monocle3::cluster_cells(cds)
  SummarizedExperiment::colData(cds)$partition <- rep(partition, dim(SummarizedExperiment::colData(cds))[1])
  #return(cds)

  cds <- monocle3::learn_graph(cds, use_partition = F, close_loop = F)

  return(cds)

}

#' Import trajectories and pseudotime from CDS to a Flow Object
#'
#'
#'
#' @param flow_object A Flow Object
#' @param cds a 'cell_data_set' object to import.
#' @param trajectory_name Name of the trajectory to add to the flow object.
#'
#' @return A Flow Object with pseudotime values for given trajectory
#'
#' @export
#'

import_trajectory <- function(flow_object,
                              cds,
                              trajectory_name){
  coll <- checkmate::makeAssertCollection()

  if(methods::is(flow_object) != "flow_object"){
    coll$push("flow_object not recognized. Please supply a valid flow object.")
  }
  if(methods::is(cds) != "cell_data_set"){
    coll$push("Error: cds value is not a valid 'cell_data_set' object from monocle3. ")
  }
  checkmate::assertCharacter(trajectory_name, any.missing = F, null.ok = F, len = 1, .var.name = "trajectory_name", add = coll)
  checkmate::reportAssertions(coll)

  pseudo_vals <- tryCatch(monocle3::pseudotime(cds),
                          error=function(cond){
                            message("Error: cds object does not have calculated pseudotime values. Run order_cells on the cds object first. ")
                            return(NA)
                          })


  if(length(intersect(trajectory_name, colnames(get_flowSet_matrix(flow_object, add_sample_id = T, annotations = T))) > 0)){
    coll$push("Error: trajectory name cannot be identical to a marker name, a metadata variable, dimension reduction or cluster annotation.")
  }
  checkmate::reportAssertions(coll)

  partition <- unique(SummarizedExperiment::colData(cds)$partition)
  if(length(partition) > 1){
    coll$push("More than one partition detected in cds object. Can only import trajectories for individual partitions ")
  } else if(is.null(flow_object$louvain$partition)){
    coll$push("No partitions in the flow_object. First run `SetPartitions()`.")
  }



  invalid_partition <- setdiff(partition, unique(flow_object$louvain$partition))
  if(length(invalid_partition) > 0 ){
    coll$push(paste("ERROR: Partition ", sQuote(unique(SummarizedExperiment::colData(cds)$partition)), " is not found in the flow_object.", sep  =""))

  }
  checkmate::reportAssertions(coll)


  out_object <- flow_object

  ica_space_df <- t(cds@principal_graph_aux$UMAP$dp_mst) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample_name") %>%
    dplyr::select(.data$sample_name, .data$UMAP1, .data$UMAP2) %>%
    dplyr::rename(prin_graph_dim_1 = "UMAP1",
                  prin_graph_dim_2 = "UMAP2")

  dp_mst <- cds@principal_graph$UMAP

  edge_df <- igraph::as_data_frame(dp_mst) %>%
    dplyr::select(source = "from", target = "to") %>%
    dplyr::left_join(.data, ica_space_df %>%
                       dplyr::select(.data$sample_name,
                                     .data$prin_graph_dim_1,
                                     .data$prin_graph_dim_2) %>%
                       dplyr::rename(source = "sample_name",
                                     source_prin_graph_dim_1 = "prin_graph_dim_1",
                                     source_prin_graph_dim_2 = "prin_graph_dim_2"),
                     by = "source") %>%
    dplyr::left_join(.data, ica_space_df %>%
                       dplyr::select(.data$sample_name,
                                     .data$prin_graph_dim_1,
                                     .data$prin_graph_dim_2) %>%
                       dplyr::rename(target = "sample_name",
                                     target_prin_graph_dim_1 = "prin_graph_dim_1",
                                     target_prin_graph_dim_2 = "prin_graph_dim_2"),
                     by = "target") %>%
    dplyr::mutate(trajectory_name = as.character(.data$trajectory_name)) %>%
    dplyr::mutate(partition = partition)

  pseudo_df <- tibble::enframe(pseudo_vals, name = "cellID", value = trajectory_name)

  if(is.null(out_object$trajectories)){
    out_object$trajectories <- list()
  }
  if(is.null(out_object$pseudotime)){
    out_object$pseudotime <- out_object$louvain %>%
      tibble::rownames_to_column("cellID") %>%
      dplyr::select(.data$cellID) %>%
      dplyr::left_join(.data, pseudo_df, by = "cellID") %>%
      tibble::column_to_rownames("cellID")
  } else {
    tmp_pseudo_df <- out_object$pseudotime %>%
      tibble::rownames_to_column("cellID")
    tmp_pseudo_df <- tmp_pseudo_df[,c("cellID", setdiff(colnames(tmp_pseudo_df), colnames(pseudo_df)))]

    out_object$pseudotime <- tmp_pseudo_df %>%
      tibble::column_to_rownames("cellID")
  }
  out_object$trajectories[[as.character(trajectory_name)]] <- edge_df


  return(out_object)


}


#' Plot a pseudotime density/relative frequency plot
#'
#' Display a density plot for pseudotime on the X-axis, and relative density on the Y-axis. Can be grouped by metadata annotations.
#'
#' @param flow_object A Flow Object
#' @param trajectory Name of the trajectory to plot.
#' @param group.by Sample or cluster annotation to separate cells by.
#' @param adjust Bin size adjustment to compute normalization over. Defaults to 1.
#' @param y Parameter to be displayed on the Y-axis. Accepted values are 'freq' or 'density'. Defaults to 'freq'.
#'
#' @return A ggplot2 object
#'
#' @export
#'

plot_grouped_pseudotime <- function(flow_object,
                                    trajectory,
                                    group.by,
                                    adjust = 1,
                                    y = "freq"){
  coll <- checkmate::makeAssertCollection()
  if (methods::is(flow_object) != "flow_object") {
    coll$push("flow_object not recognized (nor container or panel). Please supply a valid flow object.")
  }


  checkmate::assertCharacter(trajectory, any.missing = F, null.ok = F, len = 1, .var.name = "trajectory", add = coll)
  checkmate::assertCharacter(group.by, any.missing = F, null.ok = F, len = 1, .var.name = "group.by", add = coll)
  checkmate::assertCharacter(y, any.missing = F, null.ok = F, len = 1, .var.name = "freq", add = coll)
  checkmate::assertChoice(y, choices = c("freq", "density"), .var.name = "y", add = coll)
  checkmate::assertNumeric(adjust, any.missing = F, null.ok = F, len = 1, lower = 0.1, upper = 3, .var.name = "adjust", add = coll)
  checkmate::reportAssertions(coll)

  if(is.null(flow_object$pseudotime)){
    coll$push("No pseudotime values stored in this flow_object. Please use `import_trajectory` first.")

  }
  checkmate::reportAssertions(coll)

  sample_metadata <- flowWorkspace::pData(flow_object$flowSet) %>% tibble::rownames_to_column("SampleID")
  df <- get_flowSet_matrix(flow_object, add_sample_id = T, annotations = T)


  numeric_columns <- colnames(df)[unlist(lapply(df, is.numeric), use.names = FALSE)]
  character_column <- colnames(df)[!unlist(lapply(df, is.numeric), use.names = FALSE)]

  if(!is.null(group.by)){
    if(!group.by %in% colnames(df)){
      coll$push(paste("Grouping factor ", group.by, " is not found in sample metadata"))
      checkmate::reportAssertions(coll)
    }
  }


  if(!group.by %in% colnames(df)){
    coll$push(paste("Error for dataset ", sQuote(flow_object$dataset),
                    ": Cannot group.by variable ", group.by,
                    ". Column not found in metadata. Options are ",
                    paste(sQuote(character_column), collapse = ";") , sep = "" ))
  } else if(group.by %in% numeric_columns){
    coll$push(paste("Error for dataset ", sQuote(flow_object$dataset),
                    ": Cannot group.by variable ",
                    group.by, ". Metadata column cannot be of numeric type. Options are ",
                    paste(sQuote(character_column), collapse = ";") , sep = "" ))

  }

  checkmate::reportAssertions(coll)

  group_n <- df %>%
    tidyr::gather(key = "grouping", value = "grouping_value", group.by) %>%
    dplyr::group_by(.data$grouping_value) %>%
    dplyr::summarise(total = dplyr::n()) %>%
    dplyr::ungroup()

  pseudo_df <- df %>%
    dplyr::left_join(flow_object$pseudotime %>%
                       tibble::rownames_to_column("cellID") %>%
                       tidyr::gather(key = "trajectory_name", value = "pseudotime", -.data$cellID), by = "cellID")



  invalid_traj <- setdiff(trajectory, names(flow_object$trajectories))
  if(length(invalid_traj) > 0){
    coll$push(paste("Invalid trajectory name: ", paste(sQuote(invalid_traj), collapse = "; "), sep = ""))
    checkmate::reportAssertions(coll)
  } else {

    pseudo_df <- pseudo_df[pseudo_df$trajectory_name == trajectory,]

  }

  pseudo_df <- pseudo_df %>%
    tidyr::gather(key = "grouping", value = "grouping_value", group.by) %>%
    dplyr::group_by(.data$grouping_value) %>%
    dplyr::mutate(count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(.data$pseudotime))
  if(is.factor(df[,group.by]) == TRUE){
    pseudo_df$grouping_value <- factor(pseudo_df$grouping_value, levels  = levels(df[,group.by]))
  }
  p <- ggplot(pseudo_df, aes(.data$pseudotime, after_stat(.data$count), color = .data$grouping_value )) +
    theme_bw() +
    geom_density(adjust = adjust) +
    xlab("Pseudotime") +
    ylab("Density") +
    ggtitle(paste(trajectory)) +
    guides(color=guide_legend(title=group.by))



  return(p)

}


#' Plot a pseudotime boxplot per group
#'
#' Display a boxplot for pseudotime on the Y-axis, cluster or sample annotations on the X-Axis.
#'
#' @param flow_object A Flow Object
#' @param trajectory Name of the trajectory to plot.
#' @param group.by Sample or cluster annotation to separate cells by.
#'
#' @return A ggplot2 object
#'
#' @export
#'

plot_pseudotime_boxplot <- function(flow_object,
                                    trajectory,
                                    group.by){
  coll <- checkmate::makeAssertCollection()
  if (methods::is(flow_object) != "flow_object") {
    coll$push("flow_object not recognized (nor container or panel). Please supply a valid flow object.")
  }


  checkmate::assertCharacter(trajectory, any.missing = F, null.ok = F, len = 1, .var.name = "trajectory", add = coll)
  checkmate::assertCharacter(group.by, any.missing = F, null.ok = F, len = 1, .var.name = "group.by", add = coll)

  checkmate::reportAssertions(coll)

  if(is.null(flow_object$pseudotime)){
    coll$push("No pseudotime values stored in this flow_object. Please use `import_trajectory` first.")

  }
  checkmate::reportAssertions(coll)

  sample_metadata <- flowWorkspace::pData(flow_object$flowSet) %>% tibble::rownames_to_column("SampleID")
  df <- get_flowSet_matrix(flow_object, add_sample_id = T, annotations = T)


  numeric_columns <- colnames(df)[unlist(lapply(df, is.numeric), use.names = FALSE)]
  character_column <- colnames(df)[!unlist(lapply(df, is.numeric), use.names = FALSE)]


  if(!group.by %in% colnames(df)){
    coll$push(paste("Grouping factor ", group.by, " is not found in sample metadata"))
    checkmate::reportAssertions(coll)
  }



  if(!group.by %in% colnames(df)){
    coll$push(paste("Error for dataset ", sQuote(flow_object$dataset),
                    ": Cannot group.by variable ", group.by,
                    ". Column not found in metadata. Options are ",
                    paste(sQuote(character_column), collapse = ";") , sep = "" ))
  } else if(group.by %in% numeric_columns){
    coll$push(paste("Error for dataset ", sQuote(flow_object$dataset),
                    ": Cannot group.by variable ",
                    group.by, ". Metadata column cannot be of numeric type. Options are ",
                    paste(sQuote(character_column), collapse = ";") , sep = "" ))

  }

  checkmate::reportAssertions(coll)

  pseudo_df <- df %>%
    dplyr::left_join(flow_object$pseudotime %>%
                       tibble::rownames_to_column("cellID") %>%
                       tidyr::gather(key = "trajectory_name", value = "pseudotime", -.data$cellID), by = "cellID") %>%
    dplyr::filter(!is.na(.data$pseudotime))


  invalid_traj <- setdiff(trajectory, names(flow_object$trajectories))
  if(length(invalid_traj) > 0){
    coll$push(paste("Invalid trajectory name: ", paste(sQuote(invalid_traj), collapse = "; "), sep = ""))
    checkmate::reportAssertions(coll)
  } else {

    pseudo_df <- pseudo_df[pseudo_df$trajectory_name == trajectory,]

  }

  outliers <- pseudo_df %>%
    dplyr::group_by_at(vars(group.by)) %>%
    dplyr::filter(.data$pseudotime > stats::quantile(.data$pseudotime, 0.75) + 1.5 * stats::IQR(.data$pseudotime) |
                    .data$pseudotime < stats::quantile(.data$pseudotime, 0.25) - 1.5 * stats::IQR(.data$pseudotime))
  p <- ggplot(pseudo_df, aes_string(x = paste("`", group.by, "`", sep = ""), y = "pseudotime", color = paste("`", group.by, "`", sep = ""))) +
    geom_boxplot(outlier.shape = NA) +
    theme_bw()



  return(p)

}

#' Plot a pseudotime marker heatmap
#'
#' Display the relative expression for each  marker relative to pseudoplot.
#'
#' @param flow_object A Flow Object
#' @param trajectory Name of the trajectory to plot.
#' @param method Clustering method to arrange markers. Defaults to "ward.D2".  Accepted values can be found with ?hclust.
#' @param markers List of markers to display. If set to NULL, displays all included (i.e. drivers of clustering) markers. 
#' @param show_control Logical argument to determine whether or not to show reference controls. Defaults to FALSE.
#' 
#' @return A ComplexHeatmap
#'
#' @export
#'


plot_pseudotime_markers <- function(flow_object,
                                    trajectory,
                                    method = "ward.D2",
                                    markers = NULL, 
                                    show_control = FALSE){
  coll <- checkmate::makeAssertCollection()
  if (methods::is(flow_object) != "flow_object") {
    coll$push("flow_object not recognized (nor container or panel). Please supply a valid flow object.")
  }

  checkmate::assertLogical(show_control, null.ok = F, any.missing = F, len = 1, .var.name = "show_control", add = coll  )
  checkmate::assertCharacter(trajectory, any.missing = F, null.ok = F, len = 1, .var.name = "trajectory", add = coll)
  checkmate::assertCharacter(markers, any.missing = F, null.ok = T, min.len = 1, .var.name = "markers", add = coll)
  
  checkmate::reportAssertions(coll)

  if(is.null(flow_object$pseudotime)){
    coll$push("No pseudotime values stored in this flow_object. Please use `import_trajectory` first.")

  }
  
  if(!is.null(markers)){
    valid_markers <- Get_MarkerNames(flow_object, show_channel_name = F, select = "all")
    invalid_markers <- setdiff(markers, valid_markers)
    if(length(invalid_markers) >0){
      coll$push(paste("Invalid marker names: ", paste(sQuote(invalid_markers), collapse = "; "),
              "."))
    }
  }
  
  if(show_control == T){
    if(is.null(flow_object$staining_controls)){
      coll$push("Flow object does not contain staining controls. Add them first using add_staining_controls(). ")
    } else if(methods::is(flow_object$staining_controls) != "flowSet"){
      coll$push("Staining controls not recognized.  Add them first using add_staining_controls(). ")
    } else if(length(setdiff(names(flow_object$staining_controls@frames),
                             c("Unstained", flowCore::markernames(flow_object$flowSet)))) > 0){
      coll$push("Staining controls not recognized.  Add them first using add_staining_controls(). ")
    }
  }
  # 
  # if(show_control == TRUE){
  #   ref_df <- get_reference_matrix(flow_object, add_sample_id = T)
  #   if(markers %in% ref_df$SampleID){
  #     ref_df <- ref_df[ref_df$SampleID == marker,]
  #   } else{
  #     ref_df <- ref_df[ref_df$SampleID == "Unstained",]
  #   }
  #   ref_df <- ref_df %>%
  #     dplyr::mutate(louvain = "Control")
  #   ref_df <- ref_df[,c("louvain", "SampleID", marker)]
  #   df <- rbind(df, ref_df)
  #   
  #   
  #   
  # }
  checkmate::reportAssertions(coll)

  trans.obj <- flowWorkspace::flowjo_biexp_trans(equal.space = TRUE)
  df <- get_flowSet_matrix(flow_object, add_sample_id = T)

  pseudo_df <-  flow_object$pseudotime %>%
    tibble::rownames_to_column("cellID") %>%
    tidyr::gather(key = "trajectory_name", value = "pseudotime", -.data$cellID) %>%
    dplyr::filter(!is.na(.data$pseudotime))

  if(!is.null(markers)){
    mat <- df[,c(markers, "cellID")] 
    
  } else {
    mat <- df[,c(flow_object$parameters$include_markers, "cellID")]
  }
  mat <- mat %>%
    tidyr::gather(key = "marker", value = "value", -.data$cellID) %>%
    dplyr::mutate(value = ifelse(.data$value < -100, -100, .data$value )) %>%
    dplyr::mutate(value = trans.obj$transform(.data$value)) %>%
    tidyr::spread(.data$marker, .data$value) %>%
    tibble::column_to_rownames("cellID")



  pseudo_df <-  flow_object$pseudotime %>%
    tibble::rownames_to_column("cellID") %>%
    tidyr::gather(key = "trajectory_name", value = "pseudotime", -.data$cellID) %>%
    dplyr::filter(!is.na(.data$pseudotime))

  invalid_traj <- setdiff(trajectory, names(flow_object$trajectories))
  if(length(invalid_traj) > 0){
    coll$push(paste("Invalid trajectory name: ", paste(sQuote(invalid_traj), collapse = "; "), sep = ""))
    checkmate::reportAssertions(coll)
  } else {
    pseudo_df <- pseudo_df[pseudo_df$trajectory_name == trajectory,] %>%
      dplyr::arrange(.data$pseudotime) %>%
      dplyr::left_join(flow_object$louvain %>%
                         tibble::rownames_to_column("cellID"), by = "cellID") %>%
      dplyr::mutate(louvain = as.character(.data$louvain)) %>%
      tibble::column_to_rownames("cellID")
    louv_levels <- as.character(sort(unique(as.numeric(as.character(pseudo_df$louvain)))))
    pseudo_df$louvain <- factor(pseudo_df$louvain, levels = louv_levels)

  }

  mat <-  mat[row.names(pseudo_df),] %>%
    tibble::rownames_to_column("cellID") %>%
    tidyr::gather(key = "marker", value = "value", -.data$cellID) %>%
    dplyr::group_by(.data$marker) %>%
    dplyr::mutate(value = dplyr::ntile(.data$value, n = 100)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$marker) %>%
    dplyr::mutate(scaled = scale(.data$value, scale = T, center = T)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.data$value) %>%
    tidyr::spread(.data$marker, .data$scaled) %>%
    tibble::column_to_rownames("cellID") #%>%
  mat <- t(mat)

  mat <- mat[,row.names(pseudo_df)]

  louv_cols <- scales::hue_pal()(length(unique(pseudo_df$louvain)))
  names(louv_cols) <- unique(pseudo_df$louvain)
  ComplexHeatmap::ht_opt("message" = F)

  hm <- ComplexHeatmap::pheatmap(mat,
                                 name = "Z-Score",
                                 color = viridis::inferno(100),
                                 #color = colorRampPalette(c("darkblue", "blue","white", "red", "darkred"))(100),
                                 scale = "none",
                                 annotation_col = pseudo_df[c("pseudotime", "louvain")],
                                 cluster_cols = F,
                                 cluster_rows = T,
                                 show_colnames = F,
                                 clustering_method = method,
                                 main = trajectory,
                                 annotation_colors = list("louvain" = louv_cols,
                                                          #"pseudotime" = c("lightgray", "darkgreen")),
                                                          "pseudotime" = viridis::viridis(100)),
                                 heatmap_legend_param = list(direction = "vertical"))
  return(hm)


}



#' Plot a pseudotime marker heatmap
#'
#' Display the relative expression for each  marker relative to pseudoplot.
#'
#' @param cds A Flow Object
#' @param trajectory Name of the trajectory to plot.
#' @return A ComplexHeatmap
#'
#' @export
#'


FindTrajectoryMarkers <- function(cds,
                                    trajectory){
  coll <- checkmate::makeAssertCollection()
  if(methods::is(cds) != "cell_data_set"){
    coll$push("Error: cds value is not a valid 'cell_data_set' object from monocle3. ")
  }
  checkmate::assertCharacter(trajectory, any.missing = F, null.ok = F, len = 1, .var.name = "trajectory", add = coll)
  #assertCharacter(group.by, any.missing = F, null.ok = F, len = 1, .var.name = "trajectory", add = coll)
  checkmate::reportAssertions(coll)
  

  
  
  # 
  # pseudo_df <-  flow_object$pseudotime %>%
  #   tibble::rownames_to_column("cellID") %>%
  #   tidyr::gather(key = "trajectory_name", value = "pseudotime", -.data$cellID) %>%
  #   dplyr::filter(!is.na(.data$pseudotime))
  # 
  # mat <- df[,c(flow_object$parameters$include_markers, "cellID")] %>%
  #   tidyr::gather(key = "marker", value = "value", -.data$cellID) %>%
  #   dplyr::mutate(value = ifelse(.data$value < -100, -100, .data$value )) %>%
  #   dplyr::mutate(value = trans.obj$transform(.data$value)) %>%
  #   tidyr::spread(.data$marker, .data$value) %>%
  #   tibble::column_to_rownames("cellID")
  # 
  # 
  # 
  # pseudo_df <-  flow_object$pseudotime %>%
  #   tibble::rownames_to_column("cellID") %>%
  #   tidyr::gather(key = "trajectory_name", value = "pseudotime", -.data$cellID) %>%
  #   dplyr::filter(!is.na(.data$pseudotime))
  # 
  # invalid_traj <- setdiff(trajectory, names(flow_object$trajectories))
  # if(length(invalid_traj) > 0){
  #   coll$push(paste("Invalid trajectory name: ", paste(sQuote(invalid_traj), collapse = "; "), sep = ""))
  #   checkmate::reportAssertions(coll)
  # } else {
  #   pseudo_df <- pseudo_df[pseudo_df$trajectory_name == trajectory,] %>%
  #     dplyr::arrange(.data$pseudotime) %>%
  #     dplyr::left_join(flow_object$louvain %>%
  #                        tibble::rownames_to_column("cellID"), by = "cellID") %>%
  #     dplyr::mutate(louvain = as.character(.data$louvain)) %>%
  #     tibble::column_to_rownames("cellID")
  #   louv_levels <- as.character(sort(unique(as.numeric(as.character(pseudo_df$louvain)))))
  #   pseudo_df$louvain <- factor(pseudo_df$louvain, levels = louv_levels)
  #   
  # }
  # 
  # mat <-  mat[row.names(pseudo_df),] %>%
  #   tibble::rownames_to_column("cellID") %>%
  #   tidyr::gather(key = "marker", value = "value", -.data$cellID) %>%
  #   dplyr::group_by(.data$marker) %>%
  #   dplyr::mutate(value = dplyr::ntile(.data$value, n = 100)) %>%
  #   dplyr::ungroup() %>%
  #   dplyr::group_by(.data$marker) %>%
  #   dplyr::mutate(scaled = scale(.data$value, scale = T, center = T)) %>%
  #   dplyr::ungroup() %>%
  #   dplyr::select(-.data$value) %>%
  #   tidyr::spread(.data$marker, .data$scaled) %>%
  #   tibble::column_to_rownames("cellID") #%>%
  # mat <- t(mat)
  # 
  # mat <- mat[,row.names(pseudo_df)]
  # 
  # louv_cols <- scales::hue_pal()(length(unique(pseudo_df$louvain)))
  # names(louv_cols) <- unique(pseudo_df$louvain)
  # ComplexHeatmap::ht_opt("message" = F)
  # 
  # hm <- ComplexHeatmap::pheatmap(mat,
  #                                name = "Z-Score",
  #                                color = viridis::inferno(100),
  #                                #color = colorRampPalette(c("darkblue", "blue","white", "red", "darkred"))(100),
  #                                scale = "none",
  #                                annotation_col = pseudo_df[c("pseudotime", "louvain")],
  #                                cluster_cols = F,
  #                                cluster_rows = T,
  #                                show_colnames = F,
  #                                clustering_method = method,
  #                                main = trajectory,
  #                                annotation_colors = list("louvain" = louv_cols,
  #                                                         #"pseudotime" = c("lightgray", "darkgreen")),
  #                                                         "pseudotime" = viridis::viridis(100)),
  #                                heatmap_legend_param = list(direction = "vertical"))
  # return(hm)
  
  
}




