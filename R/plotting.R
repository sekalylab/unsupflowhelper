#' @include generics.R
#' @include utilities.R
#' @import ggplot2
#' @import magrittr
#' @import ggpubr
#' @import ggridges
#' @import ggrepel
#' @import flowCore
#' @importFrom rlang .data

NULL

#'
#' Plot Cluster histogram
#'
#' Graphs a histogram, similar to a 1D FlowJo, plot of the expression of a selected marker for each cluster in a Flow Object.
#'
#' @param flow_object A Flow Object
#' @param marker Select which marker to make histograms for.
#' @param alpha Control point transparency, 1 is opaque, 0 is invisible. Defaults to 0.5
#' @param clusters Filter on a subset of louvain clusters if given a vector of clusters. NULL by default.
#' @param limits Sets the lower and upper limits of the x axis. Defaults to `c(-100, 400000)`
#' @param show_control Determine if the histogram for associated control is added to the plot. Requires prior addition of controls with `add_staining_controls()`. Defaults to FALSE
#'
#' @return A ggplot object
#'
#' @export
#'

plot_cluster_histograms <- function(flow_object,
                                    marker,
                                    clusters = NULL,
                                    limits = c(-100, 400000),
                                    alpha = 0.5,
                                    show_control = FALSE){
  coll <- checkmate::makeAssertCollection()
  if(methods::is(flow_object) != "flow_object"){
    coll$push("flow_object not recognized. Please supply a valid flow object.")
  }

  checkmate::assertNumeric(alpha, lower = 0.01, upper = 1, len = 1, add = coll , any.missing = F, null.ok = F, .var.name = "alpha")
  checkmate::assertNumeric(limits, lower = -1000, upper = 450000, len = 2, add = coll , any.missing = F, null.ok = F, .var.name = "limits")
  checkmate::assertLogical(show_control, len = 1, add = coll, any.missing = F, null.ok = F, .var.name = "show_control" )

  if(!"louvain" %in% colnames(get_flowSet_matrix(flow_object))){
    coll$push(paste("Louvain clusters not found in flow object. First run ",sQuote("Cluster_flow"), ". ", sep = ""))
  }
  checkmate::reportAssertions(coll)

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

  checkmate::assertChoice(marker, choices =  flowCore::markernames(flow_object$flowSet), add = coll, null.ok = F, .var.name = "marker" )
  checkmate::assertCharacter(marker, len = 1, add = coll, any.missing = F ,null.ok = F, .var.name = "marker" )
  checkmate::reportAssertions(coll)


  flowDF_annot <- get_flowSet_matrix(flow_object, add_sample_id = T)

  flowDF_annot <- flowDF_annot[,c("louvain", "SampleID", marker)]
  #flowDF_annot <- flowDF_annot[,intersect(c("louvain", "SampleID","cellID", marker), colnames(flowDF_annot))]
  df <- flowDF_annot
  #return(df)

  if(!is.null(clusters)){
    invalid_filter <- setdiff(as.character(clusters), as.character(unique(df$louvain)))
    if(length(invalid_filter) > 0){
      coll$push(paste("Invalid louvain cluster filter values :" ,
                      paste(invalid_filter, collapse = ";"),
                      ". Valid values range from ", sQuote(min(df$louvain)), " to ",
                      sQuote(max(df$louvain)), ". ", sep = ""))
    }
    df <- df[df$louvain %in% as.numeric(clusters),]

  }
  louv_levels <- as.character(sort(unique(as.numeric(df$louvain))))
  #return(louv_levels)
  if(show_control == TRUE){
    ref_df <- get_reference_matrix(flow_object, add_sample_id = T)
    if(marker %in% ref_df$SampleID){
      ref_df <- ref_df[ref_df$SampleID == marker,]
    } else{
      ref_df <- ref_df[ref_df$SampleID == "Unstained",]
    }
    ref_df <- ref_df %>%
      dplyr::mutate(louvain = "Control")
    ref_df <- ref_df[,c("louvain", "SampleID", marker)]
    df <- rbind(df, ref_df)


    df$louvain <- factor(df$louvain, levels = c("Control", louv_levels))
  } else{
    df$louvain <- factor(df$louvain, levels = louv_levels)
  }
  p <- ggplot(df, aes_string(x = paste("`", marker, "`", sep = ""), y = "louvain", fill = "louvain", color = "louvain")) +
    ggridges::geom_density_ridges(alpha = alpha) +
    theme_bw() +
    ggcyto::scale_x_flowjo_biexp() +
    coord_cartesian(xlim = limits)
  return(p)
}


#'
#' Plot Cluster jitter plot
#'
#' Graphs a jitter plot of the frequency of cells in each sample for each cluster/cluster annotation.
#' This can then be set to be grouped by a metadata variable
#'
#' @param flow_object A Flow Object
#' @param annotation Cluster annotation to display. Defaults to `louvain`.
#' @param value Which measure of abundance to display on the plot: can acdept either `relative` or `absolute` frequencies. Defaults to `relative`.
#' @param louvain_select Filter on a subset of louvain clusters if given a vector of clusters. `NULL`by default.
#' @param group Pass a metadata variable to compare frequencies across. Defaults to `NULL`.
#' @param stats Selection of a statistical test to compare groups (if using the `groups` parameter). Accepts `wilcox`, `t.test` or `none`. Defaults to `wilcox`.
#'
#' @return A ggplot object
#'
#' @export
#'
cluster_jitterplot <- function(flow_object,
                               annotation = "louvain",
                               value = "relative",
                               louvain_select = NULL,
                               group = NULL,
                               stats = NULL){
  coll <- checkmate::makeAssertCollection()
  if(methods::is(flow_object) != "flow_object"){
    coll$push("flow_object not recognized. Please supply a valid flow object.")
  }
  checkmate::assertCharacter(as.character(louvain_select), any.missing = F, null.ok = T ,.var.name = "louvain_select", add = coll)
  checkmate::assertCharacter(annotation, null.ok = F ,any.missing = F, len = 1, .var.name = "annotation", add = coll)
  checkmate::assertCharacter(value, null.ok = F ,any.missing = F, len = 1, .var.name = "value", add = coll)
  checkmate::assertChoice(value, choices = c("relative", "absolute"), .var.name = "value", add = coll)
  checkmate::assertCharacter(group, null.ok = T ,any.missing = F, len = 1, .var.name = "group", add = coll)
  checkmate::assertCharacter(stats, null.ok = T ,any.missing = F, len = 1, .var.name = "stats", add = coll)
  if(!is.null(stats)){
    checkmate::assertChoice(stats, choices = c("wilcox", "t.test", "anova", "kruskall"), .var.name = "stats", add = coll)
  }
  checkmate::reportAssertions(coll)

  df <- get_cluster_frequencies(flow_object,
                                louvain_select = louvain_select,
                                annotation = annotation)

  if(annotation == "louvain"){
    louv_levels <- as.character(sort(unique(as.numeric(as.character(df$louvain)))))
    df$louvain <- factor(df$louvain ,levels = louv_levels)
  }


  if(value == "relative"){
    y_col <- "relative_frequency"
    y_title <- "Relative frequency (%)"
  } else {
    y_col <- "absolute_counts"
    y_title <- "Counts"
  }

  if(!is.null(group)){
    if(!group %in% colnames(flowWorkspace::pData(flow_object$flowSet))){
      coll$push(paste("Grouping factor ", group, " is not found in sample metadata"))
      checkmate::reportAssertions(coll)
    } else {
      numeric_columns <- colnames(df)[unlist(lapply(df, is.numeric), use.names = FALSE)]
      character_column <- colnames(df)[!unlist(lapply(df, is.numeric), use.names = FALSE)]

     if(group %in% numeric_columns){
      coll$push(paste("Error for dataset ", sQuote(flow_object$dataset),
                      ": Cannot group.by variable ",
                      group, ". Metadata column cannot be of numeric type. Options are ",
                      paste(sQuote(character_column), collapse = ";") , sep = "" ))

     }
    checkmate::reportAssertions(coll)

      p <- ggplot(df, aes_string(x = group, y = y_col, color = paste("`", group, "`", sep = ""))) +
        geom_boxplot() +
        geom_jitter() +
        ylab(y_title) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust=1))
      if(!is.null(stats)){
        p <- p + ggpubr::stat_compare_means(method = stats)
      }
      p <- p + facet_wrap(stats::as.formula(paste("~ `", annotation,"`", sep ="")), scales = "free")
    }

  } else {
    p <- ggplot(df, aes_string(x = paste("`", annotation, "`", sep = ""), y = y_col, color = paste("`", annotation, "`", sep = ""))) +
      geom_boxplot() +
      geom_jitter() +
      theme_bw() +
      ylab(y_title) +
      theme(axis.text.x = element_text(angle = 45, hjust=1))
  }

  return(p)

}


#'
#' Plot Stacked barplot
#'
#' Graphs a histogram, similar to a 1D FlowJo plot, of the expression of a selected marker for each cluster in a Flow Object.
#'
#' @param flow_object A Flow Object
#' @param annotation Cluster annotation to display. Defaults to `louvain`.
#' @param group Pass a metadata variable to compare frequencies across. Defaults to `NULL`.
#'
#' @return A ggplot object
#'
#' @export
#'


cluster_stacked_barplot <- function(flow_object,
                                    annotation = "louvain",
                                    group = NULL){
  coll <- checkmate::makeAssertCollection()
  if(methods::is(flow_object) != "flow_object"){
    coll$push("flow_object not recognized. Please supply a valid flow object.")
  }
  checkmate::assertCharacter(annotation, null.ok = F ,any.missing = F, len = 1, .var.name = "annotation", add = coll)
  checkmate::assertCharacter(group, null.ok = T ,any.missing = F, len = 1, .var.name = "group", add = coll)
  checkmate::reportAssertions(coll)

  df <- get_cluster_frequencies(flow_object,
                                annotation = annotation)

  if(annotation == "louvain"){
    louv_levels <- as.character(sort(unique(as.numeric(as.character(df$louvain)))))
    df$louvain <- factor(df$louvain ,levels = louv_levels)
  }


  p <- ggplot(df, aes_string(x = "SampleID", y = "relative_frequency", fill = paste("`", annotation, "`", sep = ""))) +
    theme_bw()  +
    geom_bar(position = "stack", stat = "identity") +
    ylab("Relative frequency (%)") +
    theme(axis.text.x = element_text(angle = 45, hjust=1))

  if(!is.null(group)){
    if(!group %in% colnames(flowWorkspace::pData(flow_object$flowSet))){
      coll$push(paste("Grouping factor ", group, " is not found in sample metadata"))
    } else {
      p <- p +
        facet_wrap(stats::as.formula(paste("~ `", group, "`", sep ="")), scales = "free")
    }
  }
  return(p)

}


#'
#' Plot Cluster Heatmap
#'
#' Graphs a histogram, similar to a 1D FlowJo plot, of the expression of a selected marker for each cluster in a Flow Object.
#'
#' @param flow_object A Flow Object
#' @param annotation A cluster annotation to display.
#' @param show_control Logical argument to determine whether or not show reference control MFI to the heatmap. Defaults to FALSE.
#' @param min_mfi Minimum MFI value to be displayed on the heatmap. Any value below this will be set to this value. Useful if you have a few extremely negative events. Defaults to -1000.
#' @param clustering_method Clustering method for hierarchical clustering. Viable options can be found with ?hclust
#' @param clustering_distance Clustering distance for hierarchical clustering. Viable options can be found with ?dist
#' @param add_mfi Logical argument to determine whether or not to add the MFI representation of each cluster on the absolute scale when set to TRUE.
#' @param louvain_select Select louvain clusters to represent on the heatmap
#' @param annotation_df Provide a custom cluster annotation data.frame.
#'
#' @return A ggplot object
#'
#' @export
#'

plot_cluster_heatmap <- function(flow_object,
                                 annotation = NULL,
                                 show_control = FALSE,
                                 min_mfi = -1000,
                                 clustering_method = "ward.D2",
                                 clustering_distance = "euclidean",
                                 louvain_select = NULL,
                                 add_mfi = FALSE,
                                 annotation_df = NULL){
  coll <- checkmate::makeAssertCollection()
  if(methods::is(flow_object) != "flow_object"){
    coll$push("flow_object not recognized. Please supply a valid flow object.")
  }
  checkmate::assertLogical(show_control, len = 1, add = coll, any.missing = F, null.ok = F, .var.name = "show_control" )
  checkmate::assertCharacter(clustering_method, len = 1, add = coll, any.missing = F, null.ok = F, .var.name = "clustering_method" )
  checkmate::assertNumeric(min_mfi,lower = -1500, upper = 0,  len = 1, add = coll, any.missing = F, null.ok = F, .var.name = "min_mfi" )
  checkmate::assertCharacter(clustering_distance, len = 1, add = coll, any.missing = F, null.ok = F, .var.name = "clustering_distance" )
  checkmate::assertCharacter(annotation, len = 1, add = coll, any.missing = F, null.ok = T, .var.name = "annotation" )
  checkmate::assertChoice(clustering_method, choices = c("ward,D", "ward.D2", "single", "complete", "average",
                                                         "mcquitty", "median", "centroid"),
                          .var.name = "clustering_method", add = coll)
  checkmate::assertChoice(clustering_distance, choices = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"),
                          .var.name = "clustering_distance", add = coll)
  checkmate::assertLogical(add_mfi, len = 1, add = coll, any.missing = F, null.ok = F, .var.name = "add_mfi" )

  flowDF_annot <- get_flowSet_matrix(flow_object, add_sample_id = T)
  if(!"louvain" %in% colnames(flowDF_annot)){
    coll$push(paste("Louvain clusters not found in flow object. First run ",sQuote("Cluster_flow"), ". ", sep = ""))
  }
  checkmate::reportAssertions(coll)

  valid_annot <- setdiff(colnames(flow_object$parameters$louvain_annotations), "louvain")
  if(!is.null(annotation)){
    if(length(valid_annot) == 0){
      coll$push("flow object does not contain cluster annotations. First run 'annotate_clusters'")
      checkmate::reportAssertions(coll)
    } else if(!annotation %in% valid_annot){
      coll$push(paste("Annotation ", sQuote(annotation), " was not found in the flow object. Valid values are ", paste(sQuote(valid_annot), collapse = "; "), ". ", sep = ""))
      checkmate::reportAssertions(coll)
    }
  }


  if(is.null(annotation_df)){
    annot_df <- flow_object$parameters$louvain_annotations
  } else{
    annot_df <- annotation_df
  }
  if(!is.null(annotation)){
    annot_df <- annot_df %>%
      dplyr::select(.data$louvain, dplyr::all_of(annotation)) %>%
      dplyr::rename(annotation = 2) %>%
      dplyr::mutate(annotation_label = ifelse(is.na(.data$annotation) | .data$annotation == "", .data$louvain,
                                              paste(.data$annotation, " (", .data$louvain,")", sep = "" )))
    flowDF_annot <- flowDF_annot %>%
      dplyr::left_join(annot_df, by = "louvain")
  } else {
    annot_df <- annot_df %>%
      dplyr::mutate(annotation_label = .data$louvain)
    flowDF_annot <- flowDF_annot %>%
      dplyr::left_join(annot_df, by = "louvain")
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

  #included_markers <- setdiff(markernames(flow_object$flowSet),flow_object$parameters$exclude_markers)
  included_markers <- flow_object$parameters$include_markers
  fDF <- flowDF_annot[,intersect(c(included_markers, "annotation_label"),colnames(flowDF_annot))] %>%
    tidyr::gather(key = "marker", value = "expression", -.data$annotation_label) %>%
    dplyr::group_by(.data$annotation_label,.data$marker) %>%
    dplyr::summarise(med = stats::median(.data$expression)) %>%
    dplyr::ungroup()

  if(show_control == TRUE){
    ref_df_tmp <- get_reference_matrix(flow_object, add_sample_id = T)
    ref_df_tmp <- ref_df_tmp[,intersect(c(included_markers, "SampleID"),colnames(ref_df_tmp))]  %>%
      tidyr::gather(key = "marker", value = "expression", -.data$SampleID) #%>%


    ref_df <- do.call("rbind", lapply(included_markers, function(x){
      df <- ref_df_tmp %>%
        dplyr::filter(.data$marker == x) %>%
        dplyr::mutate(keep = ifelse(x %in% ref_df_tmp$SampleID, "FMO", "Unstained")) %>%
        dplyr::mutate(keep2 = ifelse(.data$keep == "FMO" & .data$SampleID == x, 1,
                                     ifelse(.data$keep == "Unstained" & .data$SampleID == "Unstained", 1, 0))) %>%
        dplyr::filter(.data$keep2 == 1)

    })) %>%
      dplyr::select(-.data$keep, -.data$keep2) %>%
      dplyr::group_by(.data$marker) %>%
      dplyr::summarise(med = stats::median(.data$expression, na.rm = T)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(annotation_label = "Control") %>%
      dplyr::select(.data$marker, .data$med, .data$annotation_label)

  }

  trans.obj <- flowWorkspace::flowjo_biexp_trans(equal.space = TRUE)
  mat <- fDF %>%
    dplyr::mutate(med = trans.obj$transform(.data$med)) %>%
    tidyr::spread(.data$marker, .data$med) %>%
    tibble::column_to_rownames("annotation_label")



  if(!is.null(louvain_select)){

    invalid_filter <- setdiff(as.character(louvain_select), as.character(annot_df$louvain))
    if(length(invalid_filter) > 0){
      coll$push("Invalid louvain cluster filter values")
    }
    label_select <- annot_df[as.character(annot_df$louvain) %in% as.character(louvain_select),]$annotation_label
    mat <- mat[as.character(label_select),]
  }

  ord <- stats::hclust( stats::dist(as.matrix(scale(mat, center = T, scale = T)),
                      method = clustering_distance),
                 method = clustering_method)$order
  ord_markers <- stats::hclust( stats::dist(as.matrix(scale(t(mat), center = T, scale = T)),
                              method = clustering_distance),
                         method = clustering_method )$order

  if(show_control == TRUE){
    fDF <- as.data.frame(rbind(fDF, ref_df)) %>%
      dplyr::mutate(Type = ifelse(.data$annotation_label == "Control", "Control", "Cluster"))

  } else {
    fDF <- fDF %>%
      dplyr::mutate(Type =  "Cluster")

  }


  fDF["med"][fDF["med"] < min_mfi ] <- min_mfi
  #return(mat)

  flowDF_scaled <- fDF %>%
    dplyr::select(-.data$Type) %>%
    dplyr::mutate(med = trans.obj$transform(.data$med)) %>%
    dplyr::group_by(.data$marker) %>%
    dplyr::mutate(scaled = scale(.data$med, center = T, scale = T)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.data$med)

  colorAssign <- function(valueVector, scale_limits = NULL, colors = c("darkblue", "blue", "white", "red", "darkred"), length.vector = 100){
    colorLS <- grDevices::colorRampPalette(colors = colors)(length.vector)
    if(is.null(scale_limits)){
      breaks <- seq(-max(abs(valueVector)), max(abs(valueVector)), length.out = length.vector)
    } else {
      breaks <- seq(scale_limits[1], scale_limits[2], length.out = length.vector)
    }

    minBreak <- which(abs(breaks - min(valueVector)) == min(abs(breaks - min(valueVector))))
    maxBreak <- which(abs(breaks - max(valueVector)) == min(abs(breaks - max(valueVector))))
    return(colorLS[minBreak:maxBreak])
  }

  if(!is.null(louvain_select)){
    fDF_unfilt <- flowDF_scaled


    flowDF_scaled <- flowDF_scaled[flowDF_scaled$annotation_label %in% c("Control", as.character(label_select)),]
    colors <- colorAssign(valueVector = flowDF_scaled$scaled,
                          scale_limits = c(min(fDF_unfilt$scaled, na.rm = T),
                                           max(fDF_unfilt$scaled, na.rm = T)))
    fDF <- fDF[fDF$annotation_label %in% c("Control", as.character(label_select)),]
  } else{

    colors <- colorAssign(valueVector = flowDF_scaled$scaled)
  }
  col_order  <- row.names(mat)[ord]


  if(show_control == TRUE){
    col_order <- c("Control",col_order)

  }


  flowDF_scaled$annotation_label <- factor(flowDF_scaled$annotation_label,
                                  levels = as.character(col_order))

  flowDF_scaled$marker <- factor(flowDF_scaled$marker,
                                 levels = colnames(mat)[ord_markers])

  fDF$marker <- factor(fDF$marker,
                       levels = colnames(mat)[ord_markers])
  fDF <- fDF %>% dplyr::left_join(annot_df, by = "annotation_label")

  p <- ggplot(flowDF_scaled, aes_string(x = "annotation_label", y = "marker", fill = "scaled" )) +
    geom_tile(color = "black",
              lwd = 0.5,
              linetype = 1) +
    theme_classic() +
    scale_fill_gradientn(name= "Z-scale normalized MFI",
                         colors = colors,
                         limits = c(-ceiling(max(abs(flowDF_scaled$scaled),
                                                 na.rm = T)),
                                    ceiling(max(abs(flowDF_scaled$scaled),
                                                na.rm = T)))) +
    xlab("Louvain cluster") +
    ylab("Markers") +
    guides(fill = guide_colourbar(barwidth = 15,
                                  barheight = 0.5)) +
    theme(legend.position="bottom",
          axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 0.95))
  if(add_mfi == T){

    p2 <- ggplot(fDF, aes_string(x = "marker", y = "med", color = "annotation_label", shape = "Type")) +
      geom_point(size = 2) +
      theme_classic() +
      ggcyto::scale_y_flowjo_biexp() +
      theme(axis.title.y = element_blank(),
            axis.text.y = element_blank()) +
      guides(color = guide_legend(ncol = 2)) +
      ylab("Cluster MFI") +
      coord_cartesian(ylim = c(-1000, 400000)) +
      coord_flip()

    return(p + p2)

  } else {
    return(p)

  }

}


