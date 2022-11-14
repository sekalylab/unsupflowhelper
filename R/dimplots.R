#' @include generics.R
#' @include utilities.R
#' @import ggplot2
#' @import magrittr
#' @import ggrepel
#' @import flowCore
#' @importFrom rlang .data

NULL
#'
#' UMAP Dimensional reduction plot
#'
#' Graphs the output of a UMAP dimension reduction on a 2D scatter plot where each point is an event,
#' positioned based on cell embeddings determined by the `RunUMAP()` function. By default, cells are not colored,
#' but can be mapped to louvain clusters, sample or cluster annotations using the group.by parameter.
#'
#' @param flow_object A Flow Object
#' @param dot_size Control size of points. higher numbers for larger dots. Defaults to 0.5.
#' @param alpha Control point transparency, 1 is opaque, 0 is invisible. Defaults to 1
#' @param group.by Color points grouping variable. Can include 1. louvain labels; 2. cluster annotations if previously provided; 3. metadata
#' @param split.by Split the plot by sample/cluster annotations, similar to group.by. Can split by up to 3 variables.
#' @param ncols Controls the number of columns of plot if splitting by one variable
#' @param label_filter Subset louvain label filter. Useful if you only wish to show a subset of labels.
#' @param label Logical argument to determine whether or not to activate labelling on the plot. FALSE turns on legend. Defaults to TRUE.
#' @param show_trajectories Logical argument to determine whether or not to display computed trajectories on the plot. Defaults to FALSE.
#'
#' @return A ggplot object
#'
#' @seealso \code{\link{plot_density_umap}} \code{\link{plot_umap_heat}} \code{\link{plot_umap_pseudotime}}
#'
#' @export
#'



plot_umap <- function(flow_object,
                      dot_size = 0.5,
                      alpha = 1,
                      group.by = NULL,
                      split.by = NULL,
                      ncols = NULL,
                      label_filter = NULL,
                      label = TRUE,
                      show_trajectories = FALSE
){
  coll <- checkmate::makeAssertCollection()
  if(methods::is(flow_object) != "flow_object") {
    coll$push("flow_object not recognized (nor container or panel). Please supply a valid flow object.")
  }


  checkmate::assertNumeric(dot_size, lower = 0.001, upper = 2, len = 1, add = coll, any.missing = F, null.ok = F, .var.name = "dot_size")
  checkmate::assertNumeric(alpha, lower = 0.01, upper = 1, len = 1, add = coll , any.missing = F, null.ok = F, .var.name = "alpha")
  checkmate::assertCharacter(split.by, null.ok = T, max.len = 3, unique = T, .var.name = "split.by", add = coll)
  checkmate::assertCharacter(group.by, null.ok = T, max.len = 1, .var.name = "group.by", add = coll)
  checkmate::assertCharacter(as.character(label_filter), null.ok = T ,.var.name = "label_filter", add = coll)
  checkmate::assertNumeric(ncols, lower = 1, len = 1, add = coll , any.missing = F, null.ok = T, .var.name = "ncols")
  checkmate::assertLogical(label, null.ok = F, unique = T, len = 1, .var.name = "label", add = coll)
  checkmate::assertLogical(show_trajectories, null.ok = F, unique = T, len = 1, .var.name = "show_trajectories", add = coll)
  checkmate::reportAssertions(coll)

  #flowDF_annot <- get_flowSet_matrix(flow_object, add_sample_id = T, annotations = T)
  df <- get_flowSet_matrix(flow_object, add_sample_id = T, annotations = T)


  if("louvain" %in% colnames(df)){
    df$louvain <- factor(as.character(df$louvain))
    label_filt <- levels(df$louvain)
  }


  if(length(intersect(c("UMAP1", "UMAP2"), colnames(df))) < 2){
    coll$push(paste("UMAP model not found in flow object. First run ",sQuote("UMAP_flow"), ". ", sep = ""))
  }


  if(!is.null(group.by)){
    if(group.by == "louvain"){
      if(length(intersect(c("louvain"), colnames(df))) == 0){
        coll$push(paste("Louvain clusters not found in flow object. First run ",sQuote("Cluster_flow"), ". ", sep = ""))
        checkmate::reportAssertions(coll)
      } else if(!is.null(label_filter)){
        invalid_clusters <- setdiff(as.character(label_filter), levels(df$louvain))
        if(length(invalid_clusters) > 0){
          coll$push(paste("Label_filters ", paste(invalid_clusters, collapse = "; "),
                          "are invalid. Valid clusters range from ",
                          as.character(min(df$louvain)), " to ",
                          as.character(max(df$louvain)),
                          # as.character(min(flowDF_annot$louvain)), " to ",
                          # as.character(max(flowDF_annot$louvain)),
                          sep = ""))

        }
        label_filt <- as.character(label_filter)
      }

    } else if(!is.null(label_filter)){
      print("`label_filter` passed without `group.by` set to `louvain`. Ignoring...")
    }
  }



  checkmate::reportAssertions(coll)

  if(!is.null(split.by)){
    if("louvain" %in% split.by){
      if(length(intersect(c("louvain"), colnames(df))) == 0){
        coll$push(paste("Louvain clusters not found in flow object. First run ",sQuote("Cluster_flow"), ". ", sep = ""))
      }
    }
  }

  if(!is.null(flow_object$parameters$louvain_annotations)){
    clust_annot <- flow_object$parameters$louvain_annotations
    if(dim(clust_annot)[2] > 1){
      df <- df %>%
        dplyr::mutate(louvain = as.character(.data$louvain)) %>%
        dplyr::left_join(clust_annot, by = "louvain") %>%
        dplyr::ungroup()

    }
  }


  if(show_trajectories == TRUE){
    if(is.null(flow_object$trajectories)){
      coll$push("Error: No trajectories found in this flow_object.")
    }
  }

  checkmate::reportAssertions(coll)

  numeric_columns <- colnames(df)[unlist(lapply(df, is.numeric), use.names = FALSE)]
  character_column <- colnames(df)[!unlist(lapply(df, is.numeric), use.names = FALSE)]

  if(!is.null(group.by)){
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
  }
  checkmate::reportAssertions(coll)

  if(!is.null(split.by)){
    split_invalid <- setdiff(split.by,colnames(df))
    split_numeric <- intersect(split.by,numeric_columns)
    if(length(split_invalid)>0){
      coll$push(paste("Error for dataset ", sQuote(flow_object$dataset),
                      ": Cannot split.by variable(s) ", paste(split_invalid, collapse = ";"),
                      ". Column(s) not found in metadata. Options are ",
                      paste(sQuote(character_column), collapse = ";") , sep = "" ))
      checkmate::reportAssertions(coll)
    } else if(length(split_numeric) > 0){
      coll$push(paste("Error for dataset ", sQuote(flow_object$dataset),
                      ": Cannot split.by variable(s) ",
                      paste(split_numeric, collapse = "; "), ". Metadata column cannot be of numeric type. Options are ",
                      paste(sQuote(character_column), collapse = ";") , sep = "" ))
      checkmate::reportAssertions(coll)
    }
  }

  p <- ggplot(as.data.frame(df), aes_string(x= "UMAP1", y= "UMAP2" )) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    guides(colour = guide_legend(override.aes = list(size=2)))


  if(!is.null(group.by)){
    p <- p +
      geom_point(size= dot_size, alpha = alpha , shape = 16, stroke = 0,
                 aes_string(color = paste("`", group.by, "`", sep = "")))
    if(show_trajectories == TRUE){
      part_valid <- names(table(get_flowSet_matrix(flow_object, add_sample_id =  T)$partition))
      traj_df <- dplyr::bind_rows(flow_object$trajectories) %>%
        dplyr::filter(.data$partition %in% part_valid)

      p <-  p + geom_segment(data = traj_df, aes_string(x = "source_prin_graph_dim_1",
                                                 y = "source_prin_graph_dim_2",
                                                 xend = "target_prin_graph_dim_1",
                                                 yend = "target_prin_graph_dim_2"),
                             color = "black")

    }

    if(!is.null(flow_object$parameters$louvain_annotations)){
      if(group.by %in% colnames(flow_object$parameters$louvain_annotations) & label == T){
        group_ind <- which(colnames(df) == group.by)
        lc.cent <- df %>%
          dplyr::filter(.data$louvain %in% label_filt) %>%
          dplyr::group_by_at(group_ind) %>%
          dplyr::select(.data$UMAP1, .data$UMAP2) %>%
          dplyr::summarize_all("median") %>%
          dplyr::ungroup()

        p <- p +
          ggrepel::geom_label_repel(aes_string(label = paste("`", group.by, "`", sep =""),
                                      color = paste("`", group.by, "`", sep = "")), data = lc.cent) +
          guides(colour = "none")
      }

    }

  } else {
    p <- p +
      geom_point(size= dot_size, alpha = alpha , shape = 16, stroke = 0)
  }


  if(!is.null(split.by)){
    if(length(split.by) == 1){
      p <- p + facet_wrap(stats::as.formula(paste("~ `", split.by, "`", sep = "")), ncol = ncols, scales = "fixed" )
    } else {
      split_2p <- paste(unlist(lapply(split.by[c(2:length(split.by))], function(x){
        return(paste("`", x, "`", sep = ""))
      })), collapse = " + ")
      p <- p + facet_grid(stats::as.formula(paste("`", split.by[1], "`", " ~ ", split_2p, sep = "")), scales = "fixed")
    }
  }

  return(p)
}


#' UMAP Marker Fluorescence Plot
#'
#' Graphs the output of a UMAP dimension reduction on a 2D scatter plot where each point is an event,
#' positioned based on cell embeddings determined by the `RunUMAP()` function. By default, cells are not colored,
#' but can be mapped to louvain clusters, sample or cluster annotations using the group.by parameter.
#'
#' @param flow_object A Flow Object
#' @param colorGradient A vector of colors to use as a gradient for expression, from low to high. Defaults to a rainbow palette.
#' @param markers Marker(s) to use to map expression values. A minimum of one marker is required.
#' @param dot_size Control size of points. higher numbers for larger dots. Defaults to 0.5.
#' @param alpha Control point transparency, 1 is opaque, 0 is invisible. Defaults to 1
#' @param ncols Number of columns to facet multiple plots by. Applicable if using more than 1 marker.
#' @param nrows Number of rows to facet multiple plots by. Applicable if using more than 1 marker.
#' @param legend_height Sets the height of the legend. Useful if using multiple markers. Defaults to 10.
#' @param show_trajectories Logical argument to determine whether or not to display computed trajectories on the plot. Defaults to FALSE.
#' @param quantile_saturation Uses quantile based saturation to limit the effect of rare very bright events can have the the color scale as a whole. Tweak this parameter depending on the rarity of those bright events.
#'
#' @return A ggplot object
#'
#' @seealso \code{\link{plot_umap}} \code{\link{plot_density_umap}} \code{\link{plot_umap_pseudotime}}
#'
#' @export


plot_umap_heat <- function(flow_object,
                           colorGradient = c("blue4", "blue", "cyan", "green", "yellow", "orange", "red"),
                           markers,
                           dot_size = 0.5,
                           alpha = 1,
                           ncols = NULL,
                           nrows= NULL,
                           legend_height = 10,
                           show_trajectories = FALSE,
                           quantile_saturation = 1000){
  coll <- checkmate::makeAssertCollection()
  if(methods::is(flow_object) != "flow_object"){
    coll$push("flow_object not recognized. Please supply a valid flow object.")
  }

  checkmate::assertNumeric(dot_size, lower = 0.001, upper = 2, len = 1, add = coll, any.missing = F, null.ok = F, .var.name = "dot_size")
  checkmate::assertNumeric(alpha, lower = 0.01, upper = 1, len = 1, add = coll , any.missing = F, null.ok = F, .var.name = "alpha")
  checkmate::assertNumeric(ncols, len = 1, upper = 20, null.ok = T, any.missing = F, .var.name = "ncols", add = coll)
  checkmate::assertNumeric(nrows, len = 1, upper = 20, null.ok = T, any.missing = F, .var.name = "nrows", add = coll)
  checkmate::assertLogical(show_trajectories, null.ok = F, unique = T, len = 1, .var.name = "show_trajectories", add = coll)
  checkmate::assertNumeric(quantile_saturation, len = 1, lower = 100, upper = 100000, null.ok = F, any.missing = F, .var.name = "quantile_saturation", add = coll)
  checkmate::reportAssertions(coll)


  checkmate::assertCharacter(markers, min.len = 1, add = coll, any.missing = F ,null.ok = F, .var.name = "markers" )
  checkmate::assertSubset(markers, choices =  flowCore::markernames(flow_object$flowSet), add = coll, .var.name = "markers" )
  checkmate::reportAssertions(coll)


  flowDF_annot <- get_flowSet_matrix(flow_object, add_sample_id = T)
  flowDF_annot <- flowDF_annot[,intersect(c("UMAP1", "UMAP2", "louvain", "SampleID",markers), colnames(flowDF_annot))]

  sample_metadata <- flowWorkspace::pData(flow_object$flowSet) %>% tibble::rownames_to_column("SampleID")
  df <- flowDF_annot %>%
    dplyr::left_join(sample_metadata, by = "SampleID") %>%
    dplyr::ungroup()

  if(!is.null(flow_object$parameters$louvain_annotations)){
    clust_annot <- flow_object$parameters$louvain_annotations
    if(dim(clust_annot)[2] > 1){
      df <- df %>%
        dplyr::mutate(louvain = as.character(.data$louvain)) %>%
        dplyr::left_join(clust_annot, by = "louvain")

    }
  }


  if(length(intersect(c("UMAP1", "UMAP2"), colnames(flowDF_annot))) < 2){
    coll$push(paste("UMAP model not found in flow object. First run ",sQuote("UMAP_flow"), ". ", sep = ""))
  }
  checkmate::reportAssertions(coll)

  if(show_trajectories == TRUE){
    if(is.null(flow_object$trajectories)){
      coll$push("Error: No trajectories found in this flow_object.")
    }
  }

  checkmate::reportAssertions(coll)


  colors <- grDevices::colorRampPalette(colorGradient)(100)
  basic_breaks <- c(-100, 0, 1e+2, 1e+3, 1e+4, 1e+5, 262143 )

  df <- df %>%
    tidyr::gather( key = "marker", value = "expression", dplyr::all_of(markers)) %>%
    dplyr::mutate(expression = ifelse(.data$expression < -100, -100, .data$expression)) %>%
    dplyr::group_by(.data$marker) %>%
    dplyr::mutate(ntile = dplyr::ntile(.data$expression, n = quantile_saturation)) %>%
    dplyr::mutate(sat_exprs = ifelse(.data$ntile == quantile_saturation, NA, .data$expression)) %>%
    dplyr::mutate(max_exprs = max(.data$sat_exprs, na.rm = T)) %>%
    dplyr::mutate(expression = ifelse(.data$expression > .data$max_exprs, .data$max_exprs, .data$expression)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.data$sat_exprs, -.data$max_exprs, -.data$ntile)


  biexp_trans_color <- scales::trans_new(name = "flowjo_biexp_color",
                                         transform = flowWorkspace::flowjo_biexp(),
                                         inverse = flowWorkspace::flowjo_biexp(inverse = T),
                                         breaks = rev(basic_breaks),
                                         domain = c(-100, max(basic_breaks)))
  plotLS <- list()

  for(i in markers){
    df_filt <- df[df$marker == i,]

    p <- ggplot(df_filt, aes_string(x = "UMAP1", y = "UMAP2", color = "expression" )) +
      geom_point(size = dot_size, shape = 16, stroke = 0, alpha = alpha) +
      theme_bw() +
      scale_color_gradientn(name= i,
                            trans = biexp_trans_color,
                            colors = colors,
                            breaks = basic_breaks,
                            labels= scales::scientific_format(),
                            guide = guide_colorbar(barwidth = 0.5,
                                                   barheight = legend_height))
    if(show_trajectories == TRUE){
      part_valid <- names(table(get_flowSet_matrix(flow_object, add_sample_id =  T)$partition))
      traj_df <- dplyr::bind_rows(flow_object$trajectories) %>%
        dplyr::filter(.data$partition %in% part_valid)

      p <-  p + geom_segment(data = traj_df, aes_string(x = "source_prin_graph_dim_1",
                                                 y = "source_prin_graph_dim_2",
                                                 xend = "target_prin_graph_dim_1",
                                                 yend = "target_prin_graph_dim_2"),
                             color = "black")

    }
    plotLS[[i]] <- p
  }

  if(length(plotLS) > 1){
    out <- ggarrange(plotlist = plotLS,
                     ncol = ncols,
                     nrow = nrows)
  } else {
    out <- plotLS[[1]]
  }

  return(out)
}


#' UMAP Density Plot
#'
#' Graphs the output of a UMAP dimension reduction on a 2D scatter plot where each point is an event,
#' positioned based on cell embeddings determined by the `RunUMAP()` function. Cells are colored by the local density.
#' Can be split into metadata annotation to compare the relative distribution of cells across biological conditions qualitatively.
#'
#' @param flow_object A Flow Object
#' @param dot_size Control size of points. higher numbers for larger dots. Defaults to 0.5.
#' @param alpha Control point transparency, 1 is opaque, 0 is invisible. Defaults to 1
#' @param split.by Split the plot by sample/cluster annotations, similar to group.by. Can split by up to 3 variables.
#' @param ncols Controls the number of columns of plot if splitting by one variable
#' @param densityn Controls the density controller setting
#' @param show_trajectories Logical argument to determine whether or not to display computed trajectories on the plot. Defaults to FALSE.
#'
#' @return A ggplot object
#'
#' @seealso \code{\link{plot_umap}} \code{\link{plot_umap_heat}} \code{\link{plot_umap_pseudotime}}
#'
#' @export

plot_density_umap <- function(flow_object,
                              dot_size = 0.5,
                              alpha = 1,
                              split.by = NULL,
                              ncols = NULL,
                              densityn = 256,
                              show_trajectories = FALSE) {

  coll <- checkmate::makeAssertCollection()
  if (methods::is(flow_object) != "flow_object") {
    coll$push("flow_object not recognized (nor container or panel). Please supply a valid flow object.")
  }


  checkmate::assertNumeric(dot_size, lower = 0.001, upper = 3, len = 1, add = coll, any.missing = F, null.ok = F, .var.name = "dot_size")
  checkmate::assertNumeric(alpha, lower = 0.01, upper = 1, len = 1, add = coll , any.missing = F, null.ok = F, .var.name = "alpha")
  checkmate::assertVector(split.by, max.len = 3, add = coll , any.missing = F, null.ok = T, .var.name = "split.by")
  checkmate::assertNumeric(ncols, len = 1, upper = 20, null.ok = T, any.missing = F, .var.name = "ncols", add = coll)
  checkmate::assertNumeric(densityn, lower = 64, upper = 1024, len = 1, any.missing = F, null.ok = F,.var.name = "densityn", add = coll )
  checkmate::assertLogical(show_trajectories, null.ok = F, unique = T, len = 1, .var.name = "show_trajectories", add = coll)
  checkmate::reportAssertions(coll)

  flowDF_annot <- get_flowSet_matrix(flow_object, add_sample_id = T)
  flowDF_annot <- flowDF_annot[,intersect(c("UMAP1", "UMAP2", "louvain", "SampleID"), colnames(flowDF_annot))]

  sample_metadata <- flowWorkspace::pData(flow_object$flowSet) %>% tibble:: rownames_to_column("SampleID")
  df <- flowDF_annot %>%
    dplyr::left_join(sample_metadata, by = "SampleID") %>%
    dplyr::ungroup()


  if(!is.null(flow_object$parameters$louvain_annotations)){
    clust_annot <- flow_object$parameters$louvain_annotations
    if(dim(clust_annot)[2] > 1){
      df <- df %>%
        dplyr::mutate(louvain = as.character(.data$louvain)) %>%
        dplyr::left_join(clust_annot, by = "louvain")

    }

  }


  if("louvain" %in% colnames(df)){
    df$louvain <- factor(as.character(df$louvain))
    label_filt <- levels(df$louvain)
  }


  if(length(intersect(c("UMAP1", "UMAP2"), colnames(flowDF_annot))) < 2){
    coll$push(paste("UMAP model not found in flow object. First run ",sQuote("UMAP_flow"), ". ", sep = ""))
  }

  checkmate::reportAssertions(coll)

  if(!is.null(split.by)){
    if("louvain" %in% split.by){
      if(length(intersect(c("louvain"), colnames(df))) == 0){
        coll$push(paste("Louvain clusters not found in flow object. First run ",sQuote("Cluster_flow"), ". ", sep = ""))
      }
    }
  }

  checkmate::reportAssertions(coll)

  if(show_trajectories == TRUE){
    if(is.null(flow_object$trajectories)){
      coll$push("Error: No trajectories found in this flow_object.")
    }
  }

  checkmate::reportAssertions(coll)
  get_density <- function(flowDF, nbin) {
    dens <- grDevices::densCols(flowDF$UMAP1, flowDF$UMAP2, nbin = nbin,
                     colramp = grDevices::colorRampPalette(c("darkblue","darkblue", rev(grDevices::rainbow(10, end = 4/6)))))
    return(dens)
  }



  numeric_columns <- colnames(df)[unlist(lapply(df, is.numeric), use.names = FALSE)]
  character_column <- colnames(df)[!unlist(lapply(df, is.numeric), use.names = FALSE)]

  if(!is.null(split.by)){
    split_invalid <- setdiff(split.by,colnames(df))
    split_numeric <- intersect(split.by,numeric_columns)
    if(length(split_invalid)>0){
      coll$push(paste("Error for dataset ", sQuote(flow_object$dataset),
                      ": Cannot split.by variable(s) ", paste(split_invalid, collapse = ";"),
                      ". Column(s) not found in metadata. Options are ",
                      paste(sQuote(character_column), collapse = ";") , sep = "" ))

    } else if(length(split_numeric) > 0){
      coll$push(paste("Error for dataset ", sQuote(flow_object$dataset),
                      ": Cannot split.by variable(s) ",
                      paste(split_numeric, collapse = "; "), ". Metadata column cannot be of numeric type. Options are ",
                      paste(sQuote(character_column), collapse = ";") , sep = "" ))
      checkmate::reportAssertions(coll)
    }

    #df <- df[!is.na(df[split.by]),]
    if(length(split.by) > 1){
      df <- df %>% tidyr::unite("split_col", split.by, sep = "_", remove = F)
    } else {
      df$split_col <- df[[split.by]]

    }
    tmp_df <- df

    min_n <- min(as.data.frame(table(tmp_df$split_col))$Freq)

    df <- do.call("rbind",lapply(unique(tmp_df$split_col), function(x){
      filt <- tmp_df[tmp_df$split_col == x,]
      filt2 <- dplyr::slice_sample(filt, n = min_n)
      filt2$density <- get_density(filt2, densityn)
      return(filt2)
    }))
    
    p <- ggplot(as.data.frame(df), aes_string(x="UMAP1", y="UMAP2", colour = "density" )) +
      theme_bw() +
      #theme(panel.grid = element_blank()) +
      geom_point(size= dot_size, alpha = alpha , shape = 16, stroke = 0) +
      scale_color_identity()

    if(show_trajectories == TRUE){
      part_valid <- names(table(get_flowSet_matrix(flow_object, add_sample_id =  T)$partition))
      traj_df <- dplyr::bind_rows(flow_object$trajectories) %>%
                  dplyr::filter(.data$partition %in% part_valid)

      p <-  p + geom_segment(data = traj_df, aes_string(x = "source_prin_graph_dim_1",
                                                       y = "source_prin_graph_dim_2",
                                                       xend = "target_prin_graph_dim_1",
                                                       yend = "target_prin_graph_dim_2"),
                             color = "black")

    }
    
    # return(p)
    if(length(split.by) == 1){
      p <- p + facet_wrap(stats::as.formula(paste("~ `", split.by, "`", sep = "")), ncol = ncols )
    } else {
      split_2p <- paste(unlist(lapply(split.by[c(2:length(split.by))], function(x){
        return(paste("`", x, "`", sep = ""))
      })), collapse = " + ")
      p <- p + facet_grid(stats::as.formula(paste("`", split.by[1],"`" , " ~ ", split_2p, sep = "")))
    }
  

  } else {
    df$density <- get_density(df, densityn)

    p <- ggplot(as.data.frame(df), aes_string(x="UMAP1", y= "UMAP2" , colour = "density")) +
      theme_bw() +
      #theme(panel.grid = element_blank()) +
      geom_point(size= dot_size, alpha = alpha) +
      scale_color_identity()

    if(show_trajectories == TRUE){
      part_valid <- names(table(get_flowSet_matrix(flow_object, add_sample_id =  T)$partition))
      traj_df <- dplyr::bind_rows(flow_object$trajectories) %>%
                  dplyr::filter(.data$partition %in% part_valid)

      p <-  p + geom_segment(data = traj_df, aes_string(x = "source_prin_graph_dim_1",
                                                         y = "source_prin_graph_dim_2",
                                                         xend = "target_prin_graph_dim_1",
                                                         yend = "target_prin_graph_dim_2"),
                             color = "black")

    }
  }
  return(p)
}



#' UMAP Pseudotime Plot
#'
#' Graphs the output of a UMAP dimension reduction on a 2D scatter plot where each point is an event,
#' positioned based on cell embeddings determined by the `RunUMAP()` function. Cells are colored by the pseudotime value along the chosen trajectory.
#' Requires a trajectory analysis to be first performed on the flow object.
#'
#' @param flow_object A Flow Object
#' @param dot_size Control size of points. higher numbers for larger dots. Defaults to 0.5.
#' @param alpha Control point transparency, 1 is opaque, 0 is invisible. Defaults to 1
#' @param trajectory Name of the trajectory from which to display pseudotime.
#' @param viridis_palette Name of the viridis palette to use to color pseudotime values.
#'
#' @return A ggplot object
#'
#' @seealso \code{\link{plot_umap}} \code{\link{plot_umap_heat}} \code{\link{plot_density_umap}}
#'
#' @export

plot_umap_pseudotime <- function(flow_object,
                                 dot_size = 0.5,
                                 alpha = 1,
                                 trajectory,
                                 viridis_palette = "viridis"){
  coll <- checkmate::makeAssertCollection()
  if (methods::is(flow_object) != "flow_object") {
    coll$push("flow_object not recognized (nor container or panel). Please supply a valid flow object.")
  }


  checkmate::assertNumeric(dot_size, lower = 0.001, upper = 2, len = 1, add = coll, any.missing = F, null.ok = F, .var.name = "dot_size")
  checkmate::assertNumeric(alpha, lower = 0.01, upper = 1, len = 1, add = coll , any.missing = F, null.ok = F, .var.name = "alpha")
  checkmate::assertCharacter(trajectory, any.missing = F, null.ok = F, len = 1, .var.name = "trajectory", add = coll)
  checkmate::assertCharacter(viridis_palette, any.missing = F, null.ok = F, len = 1, .var.name = "viridis_palette", add = coll)
  checkmate::assertChoice(viridis_palette, choices = c("viridis", "magma", "inferno", "plasma", "cividis",
                                                       "rocket", "mako", "turbo"), .var.name = "viridis_palette", add = coll)
  checkmate::reportAssertions(coll)

  if(is.null(flow_object$pseudotime)){
    coll$push("No pseudotime values stored in this flow_object. Please use `import_trajectory` first.")

  }
  checkmate::reportAssertions(coll)





  pseudo_df <- flow_object$pseudotime %>%
    tibble:: rownames_to_column("cellID") %>%
    tidyr::gather(key = "trajectory_name", value = "pseudotime", -.data$cellID) %>%
    dplyr::left_join(flow_object$dims %>% tibble::rownames_to_column("cellID"), by = "cellID")




  invalid_traj <- setdiff(trajectory, names(flow_object$trajectories))
  if(length(invalid_traj) > 0){
    coll$push(paste("Invalid trajectory name: ", paste(sQuote(invalid_traj), collapse = "; "), sep = ""))
    checkmate::reportAssertions(coll)
  } else {
    pseudo_df <- pseudo_df[pseudo_df$trajectory_name == trajectory,]
  }


  traj_df <- flow_object$trajectories[[trajectory]]



  p <- ggplot(pseudo_df, aes_string(x="UMAP1", y="UMAP2", color = "pseudotime" )) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_point(size= dot_size , shape = 16, stroke = 0) +
    geom_segment(data = traj_df, aes_string(x = "source_prin_graph_dim_1",
                                          y = "source_prin_graph_dim_2",
                                         xend = "target_prin_graph_dim_1",
                                         yend = "target_prin_graph_dim_2"),
                 color = "black") +
    viridis::scale_color_viridis(option = viridis_palette, alpha = alpha) +
    ggtitle(trajectory)


  return(p)

}

