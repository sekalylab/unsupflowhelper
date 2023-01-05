#' @include generics.R
#' @import flowCore
#' @import magrittr
#' @import shiny
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
                               scale = FALSE){
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
  #return(cds)
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
                          scale = TRUE){
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
      dplyr::left_join(part_df, by = "cellID") %>%
      tibble::column_to_rownames("cellID")


  } else{
    part_df <- as.data.frame(partitions) %>%
      tibble::rownames_to_column("cellID") %>%
      dplyr::rename(partition = 2)

    out_object$louvain <- flow_object$louvain %>%
      tibble::rownames_to_column("cellID") %>%
      dplyr::select(.data$cellID, .data$louvain) %>%
      dplyr::left_join(part_df, by = "cellID") %>%
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
#' @param scale Logical argument to determine whether or not to scale the data prior to automatic detection of partitions.
#'
#' @return A Flow Object with partition annotations.
#'
#' @export
#'

flow_to_cds <- function(flow_object,
                        scale = FALSE){

  coll <- checkmate::makeAssertCollection()

  if(methods::is(flow_object) != "flow_object"){
    coll$push("flow_object not recognized. Please supply a valid flow object.")
  }
  #checkmate::assertCharacter(partition, any.missing = F, null.ok = F, len = 1, .var.name = "partition", add = coll)
  checkmate::assertLogical(scale, any.missing = F, null.ok = F, len = 1, .var.name = "scale", add = coll)
  #checkmate::assertNumeric(downsample, any.missing = F, null.ok = F, len = 1,lower = 500, upper = 100000, .var.name = "downsample", add = coll)
  #checkmate::assertNumeric(min, any.missing = F, null.ok = F, len = 1,lower = 100, upper = 100000, .var.name = "min", add = coll)

  checkmate::reportAssertions(coll)
  # if(is.null(flow_object$louvain$partition)){
  #   coll$push("No partitions detected in the flow_object. Run `SetPartitions()` first. ")
  # }
  # checkmate::reportAssertions(coll)
  # if(!is.null(partition)){
  # invalid.parts <- setdiff(as.character(partition), as.character(unique(flow_object$louvain$partition)))
  # if(length(invalid.parts) > 0){
  #   coll$push(paste("Invalid partition: ", paste(sQuote(invalid.parts), collapse = "; "),". Valid values are ",
  #                   paste(sQuote(unique(flow_object$louvain$partition)), collapse= "; ")))
  # }
  # checkmate::reportAssertions(coll)
  # parts <- partition
  # } else {
  #   parts <- unique(flow_object$louvain$partition)
  # }
  #return(parts)
  #out_object <- flow_object

  # if(is.null(out_object$trajectories)){
  #   out_object$trajectories <- list()
  #   out_object$louvain$pseudotime <- rep(NA, dim(get_flowSet_matrix(flow_object, add_sample_id = T))[1])
  # }
#
#   if(length(unique(flow_object$louvain$partition)) == 1){
#     message("Warning: Generating graph on whole dataset.. Consider using SetPartitions if the data is too disjointed.")
#     cds <- flow_object_to_cds(flow_object, scale = scale)
#     return(NA)
#   } else{
#     return(SubsetFlowObject(flow_object = flow_object, subset = partition == parts))
#     cds <- flow_object_to_cds(SubsetFlowObject(flow_object = flow_object, subset = partition %in% parts), scale = scale)
#   }
  cds <- flow_object_to_cds(flow_object, scale = scale)

  cds <- monocle3::cluster_cells(cds)


  clusterLS <- flow_object$louvain$louvain

  names(clusterLS) <- row.names(flow_object$louvain)
  clusterLS <- clusterLS[names(cds@clusters@listData$UMAP$clusters)]
  louv_levels <- as.character(sort(as.numeric(unique(unlist(clusterLS, use.names = F)))))
  clusterLS <- factor(clusterLS,  levels = louv_levels)
  cds@clusters@listData$UMAP$clusters <- clusterLS
  part_df <- flow_object$louvain[row.names(SummarizedExperiment::colData(cds)),]$partition
  #return(cds)
  SummarizedExperiment::colData(cds)$partition <- part_df
  return(cds)

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
#' @param
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
  if(!"cell_data_set" %in% methods::is(cds) ){
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

  # partition <- unique(SummarizedExperiment::colData(cds)$partition)
  # if(length(partition) > 1){
  #   coll$push("More than one partition detected in cds object. Can only import trajectories for individual partitions ")
  # } else if(is.null(flow_object$louvain$partition)){
  #   coll$push("No partitions in the flow_object. First run `SetPartitions()`.")
  # }
  #

#
#   invalid_partition <- setdiff(partition, unique(flow_object$louvain$partition))
#   if(length(invalid_partition) > 0 ){
#     coll$push(paste("ERROR: Partition ", sQuote(unique(SummarizedExperiment::colData(cds)$partition)), " is not found in the flow_object.", sep  =""))
#
#   }
#   checkmate::reportAssertions(coll)


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
    dplyr::left_join(ica_space_df %>%
                       dplyr::select(.data$sample_name,
                                     .data$prin_graph_dim_1,
                                     .data$prin_graph_dim_2) %>%
                       dplyr::rename(source = "sample_name",
                                     source_prin_graph_dim_1 = "prin_graph_dim_1",
                                     source_prin_graph_dim_2 = "prin_graph_dim_2"),
                     by = "source") %>%
    dplyr::left_join(ica_space_df %>%
                       dplyr::select(.data$sample_name,
                                     .data$prin_graph_dim_1,
                                     .data$prin_graph_dim_2) %>%
                       dplyr::rename(target = "sample_name",
                                     target_prin_graph_dim_1 = "prin_graph_dim_1",
                                     target_prin_graph_dim_2 = "prin_graph_dim_2"),
                     by = "target") %>%
    dplyr::mutate(trajectory_name = as.character(trajectory_name)) #%>%
   # dplyr::mutate(partition = partition)

  pseudo_df <- tibble::enframe(pseudo_vals, name = "cellID", value = trajectory_name)
  #return(pseudo_df)

  if(is.null(out_object$trajectories)){
    out_object$trajectories <- list()
  }
  if(is.null(out_object$pseudotime)){
    out_object$pseudotime <- out_object$louvain %>%
                          tibble::rownames_to_column("cellID") %>%
                          dplyr::select(.data$cellID) %>%
                          dplyr::left_join(pseudo_df, by = "cellID") %>%
                          tibble::column_to_rownames("cellID")
  } else {
    tmp_pseudo_df <- as.data.frame(out_object$pseudotime) %>%
                       tibble::rownames_to_column("cellID")

    tmp_pseudo_df <- tmp_pseudo_df[,c("cellID", setdiff(colnames(tmp_pseudo_df), colnames(pseudo_df)))] %>%
                        dplyr::left_join(pseudo_df, by = "cellID")

    out_object$pseudotime <- as.data.frame(tmp_pseudo_df) %>%
                      tibble::column_to_rownames("cellID")
  }
  out_object$trajectories[[as.character(trajectory_name)]] <- edge_df


  return(out_object)


}

#' Select an alpha value for an alpha convex hull for a given trajectory
#'
#' @param plot A Pseudotime UMAP plot. generated by the `impute_pseudotime() function`
#'
#' @return An alpha value
#'

alpha_hull_find <- function(plot){

  pseudoPlot <- plot
  pseudoDF <- pseudoPlot$data


  pseudoDF_filt <- pseudoDF %>%
    dplyr::filter(!is.na(.data$pseudotime)) %>%
    dplyr::select(.data$cellID, .data$UMAP1, .data$UMAP2) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames("cellID")



  ui <- fluidPage(
    headerPanel(
      "Impute pseudotime value within a trajectory"
    ),
    sidebarLayout(
      sidebarPanel(
      tags$h5("Select an alpha value to capture the area of values"),
       # custom column name
       sliderInput(inputId = "alpha", label = "Alpha Hull value", min = 0.02, max = 0.3, step = 0.02, value = 0.16),
       hr(),
       actionButton(inputId = "impute", label = "Impute pseudotime")
      ),
      mainPanel(
             # title for page
             # datatable output
             tags$h4("Hull area"),
             plotOutput("umap", height = "800px", width = "900px")
             # add button to finish adding column and variables
             #DT::dataTableOutput("cluster_df"),

      )
    )
  )



  server <- function(input, output, session){

    rv <- reactiveValues(ashape = alphahull::ashape(x = pseudoDF_filt,
                                                    alpha = 0.16)$edges)

    observeEvent(input$alpha,{
      rv[["ashape"]] <- as.data.frame(alphahull::ashape(x = pseudoDF_filt,
                                  alpha = input$alpha)$edges)
    })

    output$umap <- renderPlot({ plot +
                  geom_segment(data = as.data.frame(rv[["ashape"]]), color = "blue",
                                             aes(x = .data$x1, xend = .data$x2, y = .data$y1, yend = .data$y2), lwd = 2)})

    observeEvent(input$impute, {
      #stopApp(as.data.frame(rv[["ashape"]]))
      stopApp(input$alpha)
    })
  }

  sel <- shiny::runApp(shinyApp(ui = ui, server = server))


}

#' Impute pseudotime values for cells within a trajectory.
#'
#' Cluster based downsampling facilitates computation of trajectories, with an obvious caveat: it only yield pseudotime values
#' on events that were kept by downsampling. If you wish to then compare pseudotime values across metadata groups, frequencies are no
#' longer reflecting true frequencies. This function first lets you isolate the location of a given trajectory using an alpha convex hull
#' algorithm by selecting an alpha value that best captures all events in the trajectory, then imputes pseudotime on missing values
#' within that defined area by getting the mean of k neighbors.
#'
#' @param flow_object A Flow Object
#' @param trajectory Name of the trajectory to add to the flow object.
#' @param k Number of neighbors to impute pseudotime values from. Defaults to 60.
#' @param scale Logical argument to determine whether or not to scale data to determine nearest neighbours. Recommended to use the same setting used in `UMAP_flow()`. Defaults to TRUE.
#'
#' @return A Flow Object with imputed pseudotime values for given trajectory
#'
#' @export
#'

impute_pseudotime <- function(flow_object,
                               trajectory,
                               k = 60,
                               scale = TRUE){

  pseudoPlot <- plot_umap_pseudotime(flow_object, trajectory = trajectory, dot_size = 1)
  pseudoDF <- pseudoPlot$data

  pseudoDF_test <- pseudoDF %>%
                    dplyr::select(.data$cellID, .data$UMAP1, .data$UMAP2) %>%
                    tibble::remove_rownames() %>%
                    tibble::column_to_rownames("cellID")
  pseudoDF_filt <- pseudoDF %>%
                    dplyr::filter(!is.na(.data$pseudotime)) %>%
                    dplyr::select(.data$cellID, .data$UMAP1, .data$UMAP2) %>%
                    tibble::remove_rownames() %>%
                    tibble::column_to_rownames("cellID")

  alpha <- alpha_hull_find(pseudoPlot)
  # return(alpha)
  # ashape <- alphahull::ashape(x = pseudoDF_filt,
  #                             alpha = alpha)$edges
  # return(ashape)
  ahull <- alphahull::ahull(x = pseudoDF_filt,
                            alpha = alpha)


  valid <- alphahull::inahull(ahull, as.matrix(pseudoDF_test))
  keep <- pseudoDF_test[valid,]

  prep <- as.data.frame(pseudoDF) %>%
          dplyr::select(-.data$trajectory_name) %>%
          tibble::remove_rownames() %>%
          tibble::column_to_rownames("cellID")

  prep <- prep[row.names(keep),]
  filt_ids <- which(row.names(prep) %in% row.names(pseudoDF_filt))

  flowDF <- transform_flow_data(flow_object, transform = flow_object$parameters$transformation)

  if(scale == TRUE){
    flowDF <- scale(flowDF, scale = T, center = T)
  }

  #print("Performing dimension reduction using UMAP...")
  umap_nn <- suppressMessages(uwot::umap(flowDF, n_neighbors = k, n_components = 2, ret_nn = T)) #%>%

  imputed <- do.call("rbind", pbapply::pblapply(row.names(keep), function(x){
     if(x %in% row.names(pseudoDF_filt)){
       pst <- prep[x,]$pseudotime
     } else{
       nn_idx <- umap_nn$nn$euclidean$idx[x,]
       nn_names <-  row.names(umap_nn$nn$euclidean$idx[nn_idx,])
       pst <- mean(prep[nn_names,]$pseudotime, na.rm = T)
     }

     out_df <- data.frame("cellID" = x, "tmp_pst" = pst)
     return(out_df)
  }))

  out_object <- flow_object

  pseudo_df_obj <- flow_object$pseudotime %>%
                    tibble::rownames_to_column("cellID") %>%
                    dplyr::left_join(imputed, by = "cellID") %>%
                    tibble::column_to_rownames("cellID")
  pseudo_df_obj[trajectory] <- pseudo_df_obj$tmp_pst
  pseudo_df_obj$tmp_pst <- NULL
  out_object$pseudotime <- pseudo_df_obj
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
#' @return A ComplexHeatmap
#'
#' @export
#'


FindTrajectoryMarkers <- function(cds){
  coll <- checkmate::makeAssertCollection()
  if(!"cell_data_set" %in% methods::is(cds)){
    coll$push("Error: cds value is not a valid 'cell_data_set' object from monocle3. ")
  }

  #assertCharacter(group.by, any.missing = F, null.ok = F, len = 1, .var.name = "trajectory", add = coll)



  pseudo_vals <- tryCatch(monocle3::pseudotime(cds),
                          error=function(cond){
                            message("Error: cds object does not have calculated pseudotime values. Run order_cells on the cds object first. ")
                            return(NA)
                          })
  checkmate::reportAssertions(coll)

  mst <- monocle3::principal_graph(cds)$UMAP

  y_to_cells <-  monocle3::principal_graph_aux(cds)$UMAP$pr_graph_cell_proj_closest_vertex %>%
    as.data.frame()
  y_to_cells$cells <- rownames(y_to_cells)
  y_to_cells$Y <- y_to_cells$V1

  # Get the root vertices
  # It is the same node as above
  root <- cds@principal_graph_aux$UMAP$root_pr_nodes

  # Get the other endpoints
  endpoints <- names(which(igraph::degree(mst) == 1))
  endpoints <- endpoints[!endpoints %in% root]

  # For each endpoint
  # cellWeights <- lapply(endpoints, function(endpoint) {
  #   # We find the path between the endpoint and the root
  #   path <- igraph::shortest_paths(mst, root, endpoint)$vpath[[1]]
  #   path <- as.character(path)
  #   # We find the cells that map along that path
  #   df <- y_to_cells[y_to_cells$Y %in% path, ]
  #   df <- data.frame(weights = as.numeric(colnames(cds) %in% df$cells))
  #   colnames(df) <- endpoint
  #   return(df)
  # }) %>% do.call(what = 'cbind', args = .data) %>%
  #   as.matrix()
  #
  # return(cellWeights)
  # rownames(cellWeights) <- colnames(cds)
  pseudotime <- matrix(monocle3::pseudotime(cds))#, ncol = ncol(cellWeights),
                       #nrow = ncol(cds), byrow = FALSE)
  cellWeights <- matrix(1, nrow = ncol(cds), ncol = 1)
  #return(pseudotime)
  sce <- tradeSeq::fitGAM(counts = as.matrix(SummarizedExperiment::assay(cds)) + abs(min(SummarizedExperiment::assay(cds))),
                pseudotime = as.matrix(pseudotime),
                cellWeights = cellWeights)
  assc_test <- tradeSeq::associationTest(sce)
  start_test <- tradeSeq::startVsEndTest(sce)
  tradeseq_out <- list("associationTest" = assc_test[order(assc_test$waldStat, decreasing = T),],
                       "startVsEndTest" = start_test[order(start_test$waldStat, decreasing = T),],
                       "fitGAM" = sce)

  return(tradeseq_out)
}


#' Plot a pseudotime boxplot
#'
#' Display mean or median pseudotime values within a trajectory across metadata groups
#'
#' @param flow_object A Flow Object
#' @param trajectory Name of the trajectory to plot.
#' @param metric Summary statistic to display. Accepts `mean` or `median`.
#' @param group.by Metadata variable to compare pseudotime.
#' @param stats Selection of a statistical test to compare groups (if using the `groups` parameter). Accepts `wilcox`, `t.test` or `none`. Defaults to `wilcox`.
#' @return A ggplot boxplot
#'
#' @export
#'

plot_pseudotime_boxplot <- function(flow_object,
                                    trajectory,
                                    metric = "mean",
                                    group.by,
                                    stats = NULL){
  coll <- checkmate::makeAssertCollection()
  if (methods::is(flow_object) != "flow_object") {
    coll$push("flow_object not recognized (nor container or panel). Please supply a valid flow object.")
  }


  checkmate::assertCharacter(trajectory, any.missing = F, null.ok = F, len = 1, .var.name = "trajectory", add = coll)
  checkmate::assertCharacter(metric, any.missing = F, null.ok = F, min.len = 1, .var.name = "metric", add = coll)
  checkmate::assertCharacter(group.by, null.ok = F ,any.missing = F, len = 1, .var.name = "group.by", add = coll)
  checkmate::assertCharacter(stats, null.ok = T ,any.missing = F, len = 1, .var.name = "stats", add = coll)
  if(!is.null(stats)){
    checkmate::assertChoice(stats, choices = c("wilcox", "t.test", "anova", "kruskal.test"), .var.name = "stats", add = coll)
  }
  checkmate::reportAssertions(coll)

  if(is.null(flow_object$pseudotime)){
    coll$push("No pseudotime values stored in this flow_object. Please use `import_trajectory` first.")

  }

  if(!metric %in% c("mean", "median")){
    coll$push("Invalid metric: accepted values are `mean` or `median`")
  }

  if(!group.by %in% colnames(flowWorkspace::pData(flow_object$flowSet))){
    coll$push(paste("Grouping factor ", group.by, " is not found in sample metadata"))

  }
  checkmate::reportAssertions(coll)





  invalid_traj <- setdiff(trajectory, names(flow_object$trajectories))
  if(length(invalid_traj) > 0){
    coll$push(paste("Invalid trajectory name: ", paste(sQuote(invalid_traj), collapse = "; "), sep = ""))
    checkmate::reportAssertions(coll)
  }

  pseudo_df <-  flow_object$pseudotime %>%
    tibble::rownames_to_column("cellID") %>%
    tidyr::gather(key = "trajectory_name", value = "pseudotime", -.data$cellID) %>%
    dplyr::filter(.data$trajectory_name == trajectory) %>%
    tidyr::separate(.data$cellID, into = c("SampleID", "cellNumber"), sep = "__") %>%
    dplyr::group_by(.data$SampleID, .data$trajectory_name) %>%
    dplyr::summarise(mean = mean(.data$pseudotime, na.rm = T),
                     median = stats::median(.data$pseudotime, na.rm = T)) %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(Get_MetaData(flow_object) %>%
                        tibble::rownames_to_column("SampleID"), by = "SampleID")


  p <- ggplot(pseudo_df, aes_string(x = group.by, y = metric, color =  paste("`", group.by, "`", sep = ""))) +
    geom_boxplot() +
    geom_jitter() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust=1))

  if(metric == "mean"){
    p <- p +
          ylab("Mean Pseudotime")
  } else if(metric == "median"){
    p <- p +
      ylab("Median Pseudotime")
  }
  if(!is.null(stats)){
    p <- p + ggpubr::stat_compare_means(method = stats)
  }

  return(p)


}


