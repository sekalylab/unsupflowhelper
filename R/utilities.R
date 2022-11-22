#' @include generics.R
#' @import flowCore
#' @import magrittr
#' @import shiny
#' @importFrom rlang .data
#'
NULL


#' Get Sample Metadata
#'
#' Retrieves sample metadata from a Flow Object
#'
#' @param flow_object a Flow Object
#'
#' @return A data.frame
#'
#' @export
#'
Get_MetaData <- function(flow_object){
  coll = checkmate::makeAssertCollection()
  if(methods::is(flow_object) != "flow_object") { coll$push("flow_object not recognized. Please supply a valid flow panel. ")}

  checkmate::reportAssertions(coll)

  return(flowCore::pData(flow_object$flowSet))
}

#' Get Marker Names
#'
#' Retrieves marker Names from a Flow Object
#'
#' @param flow_object a Flow Object
#'
#' @return A marker name vector
#'
#' @export
#'
#'
Get_ClusterAnnotations <- function(flow_object){
    coll = checkmate::makeAssertCollection()
    if(methods::is(flow_object) != "flow_object") { coll$push("flow_object not recognized. Please supply a valid flow panel. ")}

    checkmate::reportAssertions(coll)

    return(flow_object$parameters$louvain_annotations)
}


#' Get Cluster Annotations
#'
#' Retrieves cluster annotations from a Flow Object
#'
#' @param flow_object a Flow Object
#'
#' @return A data.frame
#'
#' @export
#'
Get_MetaData <- function(flow_object){
  coll = checkmate::makeAssertCollection()
  if(methods::is(flow_object) != "flow_object") { coll$push("flow_object not recognized. Please supply a valid flow panel. ")}

  checkmate::reportAssertions(coll)

  return(flowCore::pData(flow_object$flowSet))
}

#' Get Marker Names
#'
#' Retrieves marker Names from a Flow Object
#'
#' @param flow_object a Flow Object
#' @param select Put a filter on which markers to display. Accepted values are 'all', 'included' or 'excluded'.
#' @param show_channel_name Logical argument to determine whether or not to also display associated channel names. Defaults to TRUE.
#'
#' @return A marker name vector
#'
#' @export
#'


Get_MarkerNames <- function(flow_object,
                            select = "all",
                            show_channel_name = TRUE){
  coll = checkmate::makeAssertCollection()
  if(methods::is(flow_object) != "flow_object") { coll$push("flow_object not recognized. Please supply a valid flow panel. ")}
  checkmate::assertCharacter(select, len = 1, any.missing = F, null.ok = F, .var.name = "select", add = coll)
  checkmate::assertLogical(show_channel_name, len = 1, any.missing = F, null.ok = F, .var.name = "show_channel_names", add = coll)
  checkmate::reportAssertions(coll)

  if(!select %in% c("all", "included", "excluded")){
    coll$push(paste("Invalid 'select' value ", sQuote(select),". Valid 'select' values are 'all', 'included' or 'excluded'.", sep = "")) }
  if(methods::is(flow_object) != "flow_object") { coll$push("flow_object not recognized. Please supply a valid flow panel. ")}
  checkmate::reportAssertions(coll)
  out <- flowCore::markernames(flow_object$flowSet)

  if(select == "included"){
    out <- out[out %in% flow_object$parameters$include_markers]
  } else if(select == "excluded"){
    out <- out[out %in% flow_object$parameters$exclude_markers]
  }
  if(show_channel_name == FALSE){
    out <- unname(out)
  }
  return(out)
}


#' Downsample a FlowSet to keep n events in each file
#'
#' Randomly selects n events from each file in a subset, in order to get fewer events than the total number across all files.
#'
#' @param flowSet A FlowSet
#' @param n A number of events to take per file. Values between 100 and 500000 are accepted.
#'
#' @return A FlowSet
#'
#'
#' @examples
#' library(flowCore)
#' data("GvHD")
#' Downsampling_FlowSet(flowSet = GvHD, n = 1000)
#'
#' @export
#'
Downsampling_FlowSet <- function(flowSet,
                                 n ){
  coll = checkmate::makeAssertCollection()
  checkmate::assertClass(flowSet, classes = "flowSet", .var.name = "flowSet", add = coll )
  checkmate::assertNumeric(n, len = 1, any.missing = F,lower = 100, upper = 500000, .var.name = "n", add = coll )
  checkmate::reportAssertions(coll)


  if(missing(n))
    n <- min(flowCore::fsApply(flowSet,nrow))
  fS <- flowCore::fsApply(flowSet, function(ff){
    if(nrow(ff) < n){
      mess <- paste("Sample ", ff@description$`TUBE NAME`, " (", basename(ff@description$FILENAME), ") only contains ",
                    nrow(ff), " events. Using all ", nrow(ff), " events", sep = "")
      print(mess)

      return(ff)
    } else {
      i <- sample(nrow(ff), size = n, replace=FALSE, prob = NULL)
      return(ff[i,])
    }

  })
  return(fS)
}


#' Convert a list of CSV files exported from FlowJo to a flowSet
#'
#' Randomly selects n events from each file in a subset, in order to get fewer events than the total number across all files.
#'
#' @param files a list of csv files
#'
#' @return A FlowSet
#'
#' @export
#'

csv_to_flowSet <- function(files){
  coll <- checkmate::makeAssertCollection()
  checkmate::assertCharacter(files, min.len = 2, any.missing = F, null.ok = F, .var.name = "files", add = coll)
  checkmate::reportAssertions(coll)
  checkmate::assertFileExists(files, extension = "csv", .var.name = "files", add = coll)
  checkmate::reportAssertions(coll)

  fileLS <- files
  names(fileLS) <- basename(files)
  ffLS <- sapply(names(fileLS), simplify = F, USE.NAMES = T, function(x){
    df <- readr::read_csv(file = fileLS[[x]],
                          show_col_types = F)

    param_df <- data.frame("name" = colnames(df),
                           "desc"  = colnames(df)) %>%
      dplyr::mutate(desc = ifelse(grepl("FSC-|SSC-|Time$", .data$name), NA, dplyr::desc)) %>%
      dplyr::mutate(range = 262144,
                    minRange = -111,
                    maxRange = 262144)
    row.names(param_df) <- paste("$P", row.names(param_df), sep = "")
    #return(param_df)
    ff <- methods::new("flowFrame",
              exprs = as.matrix(df),
              parameters = Biobase::AnnotatedDataFrame(param_df))
    return(ff)
  })
  return(methods::as(ffLS, "flowSet"))
}


#' Subset a FlowObject.
#'
#' Subset a FlowObject on metadata criteria, cluster annotations, dimensions or marker fluorescence.
#'
#' @param flow_object A flow_object
#' @param subset criteria for filtering. See Details.
#'
#' @return A flow_object
#'
#' @details
#' Any parameter embedded in the flow_object can be potentially used as a filter, as a logical expression.
#' This includes marker fluorescence expression, cluster name or annotation, sample metadata values, UMAP coordinates.
#' Logical expressions can be chained together using the & (AND) or | operators (OR) .
#'
#' Expression operators such as == (equals to), != (does not equal), > or < (greater or smaller than), %in% (is a part of ...) are valid.
#'
#' @export
#'


SubsetFlowObject <- function(flow_object,
                             subset){

  coll = checkmate::makeAssertCollection()
  if(methods::is(flow_object) != "flow_object") { coll$push("flow_object not recognized. Please supply a valid flow panel. ")}

  checkmate::reportAssertions(coll)

  flowDF_annot <- get_flowSet_matrix(flow_object, add_sample_id = T, annotations = T)

  expr <- if (tryCatch(expr = rlang::is_quosure(x = subset),
                       error = function(...) FALSE)) {
    expression
  }
  else if (is.call(x = rlang::enquo(arg = subset))) {
    rlang::enquo(arg = subset)
  }
  else {
    parse(text = subset)
  }

  data.subset <- flowDF_annot %>%
                dplyr::mutate(tmp = .data$cellID) %>%
                tibble::column_to_rownames("tmp")

  cells <- row.names(data.subset)[rlang::eval_tidy(expr = expr, data = data.subset)]

  filt_df <- flowDF_annot[flowDF_annot$cellID %in% cells,]
  tmp_fS <- flowSet_to_list(flow_object$flowSet)

  for(i in names(tmp_fS)){

    if(length(intersect(row.names(flowCore::exprs(tmp_fS[[i]])), filt_df$cellID)) == 0){
      tmp_fS[[i]] <- NULL
    } else{
      flowCore::exprs(tmp_fS[[i]]) <- flowCore::exprs(tmp_fS[[i]])[row.names(flowCore::exprs(tmp_fS[[i]])) %in% filt_df$cellID,]
    }

  }

  fs_new <- methods::as(tmp_fS, "flowSet")
  #return(fs_new)
  if(dim(flowWorkspace::pData(flow_object$flowSet))[2] > 1){
    pData <- flowWorkspace::pData(flow_object$flowSet)[flowWorkspace::sampleNames(fs_new),]
    #pData <- pData[row.names(pData) %in% flowWorkspace::sampleNames(fs_new),]
    flowWorkspace::pData(fs_new) <- pData
  }

  out_object <- flow_object
  out_object$flowSet <- fs_new

  if(!is.null(flow_object$dims)){
    filt_dims <- flow_object$dims[cells,]
    out_object$dims <- filt_dims
  }
  if(!is.null(flow_object$louvain)){
    louv_df <- flow_object$louvain
    if(dim(louv_df)[2] == 1){
      filt_louvain <- data.frame("louvain" = flow_object$louvain[cells,])
      row.names(filt_louvain) <- cells
    } else{
      filt_louvain <- flow_object$louvain[cells,]
    }


    out_object$louvain <- filt_louvain

  }
  if(!is.null(flow_object$pseudotime)){
    pseudo_df <- as.data.frame(flow_object$pseudotime) %>%
      tibble::rownames_to_column("cellID") %>%
      dplyr::filter(.data$cellID %in% cells) %>%
      tibble::column_to_rownames("cellID")

    out_object$pseudotime <- pseudo_df

  }

  return(out_object)

}

#' Subsample on a per cluster basis
#'
#' Subset a FlowObject on metadata criteria, cluster annotations, dimensions or marker fluorescence.
#'
#' @param flow_object A flow_object
#' @param n number of cells to subsample per clusters. Defaults to 10000
#'
#' @return A flow_object
#'
#' @export
#'


ClusterSubsampling <- function(flow_object,
                               n = 10000){
  coll = checkmate::makeAssertCollection()

  if(methods::is(flow_object) != "flow_object") { coll$push("flow_object not recognized. Please supply a valid flow panel. ")}
  checkmate::assertNumeric(n, lower = 1000, upper = 50000, len = 1, any.missing = F, null.ok = F, .var.name = "n", add = coll)
  checkmate::reportAssertions(coll)

  df <- get_flowSet_matrix(flow_object, add_sample_id = T)

  if(length(intersect(c("louvain"), colnames(df))) == 0){
    coll$push(paste("Louvain clusters not found in flow object. First run ",sQuote("Cluster_flow"), ". ", sep = ""))
    checkmate::reportAssertions(coll)
  }

  ncells <- df %>%
    dplyr::group_by(.data$louvain) %>%
    dplyr::summarise(subsample = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(subsample = ifelse(.data$subsample > n, n, .data$subsample))

  cells <- df %>%
    dplyr::group_by(.data$louvain) %>%
    tidyr::nest() %>%
    dplyr::left_join(ncells, by = "louvain" ) %>%
    dplyr::mutate(Sample = purrr::map2(.data$data, .data$subsample, dplyr::sample_n)) %>%
    tidyr::unnest(.data$Sample) #%>%


  sub <- SubsetFlowObject(flow_object, subset = .data$cellID %in% cells$cellID)

  return(sub)
}

#' Find optimal cofactors for hyperbolic sine transform
#'
#' Incrementally increase cofactor values for an hyperbolic sine transformation for each marker until the fluorescence values are not longer split into a bimodal dsitribution centered around zero.
#' This is key to avoid having clustering and dimension reduction driven by very low expression markers.
#'
#' @param flow_object A Flow Object
#'
#' @return A flow_object
#'
#' @export
#'

find_cofactors <- function(flow_object){
  coll <- checkmate::makeAssertCollection()
  if(methods::is(flow_object) != "flow_object") { coll$push("flow_object not recognized. Please supply a valid flow panel. ")}
  checkmate::reportAssertions(coll)


  flowDF <- get_flowSet_matrix(flow_object)

  #include <- setdiff(colnames(flowDF),c(flow_object$parameters$exclude_markers, "UMAP1", "UMAP2", "louvain"))
  include <- flow_object$parameters$include_markers

  flowDF <- flowDF[,include]
  markers <- pbapply::pbsapply( colnames(flowDF),simplify = F, USE.NAMES = T, function(x){
    cofactor_check <- TRUE
    cofactor <- 50
    while(cofactor_check & cofactor <= 10000){
      expr <- asinh(flowDF[,x]/cofactor)
      if(length(expr) > 0){
        quant <- stats::quantile(expr, probs = c(0.001))
        if(quant < -0.5){
          cofactor <- cofactor + 5
        } else {
          cofactor_check <- F
        }
      }
    }
    return(cofactor)
  })
  out_object <- flow_object
  out_object$parameters$cofactors <- unlist(markers, use.names = T)
  out_object$parameters$cofactor_find <- TRUE

  return(out_object)

}

#' Annotate louvain clusters
#'
#' Add annotations to louvain clusters. Multiple concurrent annotations can be used in parallel. For instance, you could have more basic cell types `(CD4 T cells, CD8 T cells, etc)`, or more granular.
#' It is possible to have multiple clusters with the same annotations if you wish to group them. Those annotations will then be useable in plots and statistics instead of the default louvain annotations.
#'
#' @param flow_object A Flow Object
#' @param cluster_annotations A data.frame containing a column with louvain labels, and additional columns with user-defined labels. Ignored if using the 'interactive' mode.
#' @param interactive Logical argument to determine whether or not to activate the interactive mode, which opens up a graphical interface to manually add labels. Defaults to FALSE.
#'
#' @return A flow_object
#'
#' @export
#'
annotate_clusters <- function(flow_object,
                              cluster_annotations = NULL,
                              interactive = FALSE){
  coll <- checkmate::makeAssertCollection()
  if (methods::is(flow_object) != "flow_object") {
    coll$push("flow_object not recognized (nor container or panel). Please supply a valid flow object.")
  }

  checkmate::assertDataFrame(cluster_annotations, null.ok = T, any.missing = F,.var.name = "cluster_annotation", add = coll)
  checkmate::assertLogical(interactive, any.missing = F, len = 1, null.ok = F, .var.name = "interactive", add = coll)
  checkmate::reportAssertions(coll)

  if(interactive == FALSE & is.null(cluster_annotations)){
    coll$push("If 'interactive' is set to FALSE, you need to provide a `cluster_annotations` data.frame. ")

  }
  if(is.null(flow_object$louvain)){
    coll$push(paste("Louvain clusters not found in flow object. First run ",sQuote("Cluster_flow"), ". ", sep = ""))
  }
  checkmate::reportAssertions(coll)

  out_object <- flow_object

  if(interactive == TRUE){
    df <- annotate_clusters_shiny(flow_object)
    out_object$parameters$louvain_annotations <- df

  } else {
    cluster_ls <- as.character(sort(unique(flow_object$louvain$louvain)))
    #if("louvain" %in% colnames(get_flowSet_matrix(flow_object))){
    if(length(intersect("louvain", colnames(cluster_annotations))) == 0){
      coll$push("Supplied `cluster_annotations` data.frame does not contain a column called `louvain`")
    } else if(length(intersect(cluster_ls, as.character(cluster_annotations$louvain))) < length(cluster_ls)){
      coll$push("Supplied `cluster_annotations` data.frame does not contain all `louvain` values in the `louvain` column")
    }

    checkmate::reportAssertions(coll)

    if(dim(flow_object$parameters$louvain_annotations)[2] == 1){
      tmp_df <- data.frame("louvain" = flow_object$parameters$louvain_annotations$louvain)
    } else {
      tmp_df <- flow_object$parameters$louvain_annotations
    }
    return(tmp_df)
    if(dim(tmp_df)[2] == 1){
      df <- cluster_annotations
    } else {
      df <- tmp_df[,c("louvain",setdiff( colnames(tmp_df), colnames(cluster_annotations) ))] %>%
        dplyr::left_join( cluster_annotations, by = "louvain")
    }

    out_object$parameters$louvain_annotations <- df
  }



  return(out_object)
}

#' Louvain annotation shiny interface
#'
#' Boots the annotate_cluster shiny interface. It is meant to be called by annotate_clusters function, not directly.
#'
#' @param flow_object A Flow Object
#'
#' @return A data.frame
#'
annotate_clusters_shiny <- function(flow_object ){
  ## adapted from https://github.com/hinkelman/Shiny-Scorekeeper/blob/master/server.R


  ref <- !is.null(flow_object$staining_controls)
  louvain_order <- levels(plot_cluster_heatmap(flow_object, show_control = ref, annotation = NULL)$data$annotation_label)

  cluster_df <- as.data.frame(flow_object$parameters$louvain_annotations)
  cluster_df2 <- as.data.frame(cluster_df[match(intersect(louvain_order, cluster_df$louvain), cluster_df$louvain),]) %>%
                  dplyr::rename(louvain = 1)
  # datatable output function
  dt_output <- function(title, id) {
    fluidRow(column(
      12,
      h4(title)), # title
      hr(),
      DT::DTOutput(id) # datatable output
    )
  }

  # render datatable function
  # render_dt = function(data, editable = list(target = "column",
  #                                            disable = list(columns = c(0))), server = TRUE, ...) {
  render_dt = function(data, editable = TRUE, server = TRUE, ...) {
    DT::renderDT(data,
                 selection = 'none',
                 server = server,
                 editable = editable,
                 rownames = F,
                 options = list(searching = F,
                                info = F,
                                iDisplayLength=50,
                                bFilter = F,
                                bPaginate = F,
                                bLengthChange = 0,
                                bInfo = 0,
                                bAutoWidth = 1,
                                ordering = F),
                 ...)
  }

  ui <- fluidPage(
    headerPanel(
      "Cluster annotations"
    ),
    fluidRow(
      column(width = 6,
             # custom column name
             textInput(inputId = "nameColumn", "New Column Name"),
             actionButton(inputId = "addColumn", "Add/Reset Column"),
             hr(),
             tags$h4("Cluster Heatmap"),
             uiOutput("column_sel"),
             plotOutput("heatmap"),
             tags$h5("Add a cluster annotation column"),
             helpText("1. Write the name of a new cluster annotation. Multiple cluster annotations can be used in parallel. ", br(),
                      "2. Click `Add Column`. ", br(),
                      "3. Add labels by double-clicking on the table cells.", br(),
                      "4. Click `Done` when you are ready to submit")


      ),

      column(width = 6,
             # title for page
             title = 'Cluster heatmap',
             # datatable output
             tags$h3("Double-click to edit table cells"),
             helpText("** Known bug: the first value entered disapears. Re-enter it and move on. "),
             dt_output('', 'x'),
             # add button to finish adding column and variables
             actionButton(inputId = "done", "Done"),
             #DT::dataTableOutput("cluster_df"),

             tags$hr()
      )
    )
  )


  server <- function(input, output, session){

    #data <- reactive({cluster_df})

    rv <- reactiveValues(data = cluster_df2)

    output$heatmap <- renderPlot({plot_cluster_heatmap(flow_object, show_control = ref,
                                                       annotation = input$column_sel,
                                                       annotation_df = rv[["data"]])})
    # observe data
    # observeEvent(data(),{
    #   reactiveData(data())
    # })


    # edit a single cell
    # make cells editable in our data
    output$x <- render_dt(rv[["data"]])

    proxy <- DT::dataTableProxy('x')

    observeEvent(input$x_cell_edit, {
      info = input$x_cell_edit
      i <- info$row
      j <- info$col + 1L  # column index offset by 1
      v <- info$value

      rv[["data"]][i, j] <- suppressWarnings(DT::coerceValue(v, rv[["data"]][i, j]))
      DT::replaceData(proxy, rv[["data"]], resetPaging = FALSE, rownames = FALSE)  # important
      # newData <- reactiveData()
      # newData[info$row, info$col] <- info$value
      # reactiveData(newData)
      #
      # replaceData(proxy, reactiveData(), resetPaging = FALSE)
    })


    # add a column
    observeEvent(input$addColumn,{
      newData <- rv[["data"]]
      newData[[input$nameColumn]] <- rep(NA, nrow(newData))
      #newData[[input$nameColumn]] <- numeric(length = nrow(newData))
      rv[["data"]] <<- newData
      DT::replaceData(proxy, rv[["data"]], resetPaging = FALSE)
      output$x = render_dt({
        rv[["data"]]
      })
    })

    observeEvent(input$validate, {
      out_df <- rv[["data"]]
      #newData <- reactiveData()
      #newData[info$row, info$col] <- info$value


    })
    output$column_sel = renderUI({
      selectizeInput("column_sel",
                     label = h5("Select Annotation Column"),
                     choices = setdiff(colnames(rv[["data"]]), "louvain"),
                     multiple = TRUE,
                     options = list(maxItems = 1))
    })

    observeEvent(input$x_cell_edit, {
      updateSelectizeInput(session, "column_sel", selected = input$select_col)
    })
    observeEvent(input$done, {
      out_df <- rv[["data"]]
      out_df <- out_df[match( cluster_df$louvain, out_df$louvain),]
      stopApp(out_df)
    })
  }

  sel <- shiny::runApp(shinyApp(ui = ui, server = server))

}


#' Add metadata annotations
#'
#' Add sample annotations to the Flow Object. This is where you add clinical variables (i.e. antibody levels, Timepoints, Treatment, etc) or technical variables (like Batch).
#' Variables can be discrete (i.e. categorical) or continuous. If using categorical, it is recommended to assign factor levels if you desire a specific order for labels to be displayed in plots, or how statistics will be performed.
#' By default, the order will be performed in alphabetical order.
#'
#' @param flow_object A Flow Object
#' @param metadata A data.frame containing sample annotations.
#' @param key_column Name of the column matching sample names.
#' @param flowobj_column Name of the column matching sample names to be matched against in the flow object.
#' @param interactive  Logical argument to determine whether or not to activate the interactive mode, which opens up a graphical interface to select a metadata file. Defaults to FALSE.
#'
#' @return A Flow Object
#'
#' @export
#'

add_sample_metadata <- function(flow_object,
                                metadata = NULL,
                                key_column = NULL,
                                flowobj_column = NULL,
                                interactive = FALSE){
  coll <- checkmate::makeAssertCollection()
  if(methods::is(flow_object) != "flow_object"){
    coll$push("flow_object not recognized. Please supply a valid flow object.")
  }
  checkmate::assertClass(metadata, classes = "data.frame", null.ok = T, .var.name = "metadata", add = coll)
  checkmate::assertCharacter(key_column, len = 1, any.missing = F, null.ok = T, .var.name = "key_column", add = coll)
  checkmate::assertCharacter(flowobj_column, len = 1, any.missing = F, null.ok = T, .var.name = "flowobj_column", add = coll)
  checkmate::reportAssertions(coll)



  checkmate::reportAssertions(coll)

  df <- as.data.frame(flowWorkspace::pData(flow_object$flowSet)) %>%
            tibble::rownames_to_column("tmp_sampleid")

  if(interactive == T){
    meta_int <- add_metadata_shiny(flow_object)
    metadata_df <- meta_int$metadata_df
    key_col <- meta_int$key
    flowobj_col <- meta_int$fS_key
  } else {
    if(is.null(metadata)){
      coll$push("Error: Must supply metadata dataframe in non-interactive mode")
      checkmate::reportAssertions(coll)
    } else if(is.null(key_column)){
      coll$push("Error: Must supply a key column in non-interactive mode")
      checkmate::reportAssertions(coll)
    } else if(is.null(flowobj_column)){
      coll$push("Error: Must supply a key flowobject column in non-interactive mode")
      checkmate::reportAssertions(coll)
    } else {
      metadata_df <- metadata
      key_col <- key_column
      flowobj_col <- flowobj_column
    }

  }

  # if(!key_col %in% colnames(metadata_df)){
  #   coll$push(paste("Sample metadata does not contain a corrresponding ", sQuote(key_col) ," column", sep = ""))
  # }
  #
  invalid_samples <- setdiff(df[[flowobj_col]], metadata_df[[key_col]])
  if(length(invalid_samples) > 0){
    coll$push(paste("Samples ", paste(sQuote(invalid_samples), collapse = "; "),
                    " were not found in the selected flow_object column.  ", sep = ""))

  }

  sample_intersect <- intersect(df[[flowobj_col]], metadata_df[[key_col]])
  if(length(sample_intersect) == 0){
    coll$push("No common samples found in the supplied metadata column vs the flow_object metadata column.")
    reportAssertions(coll)
  } else if(length(sample_intersect) < length(df[[flowobj_col]])){
    message("Warning: supplied metadata has values for only a subset of samples. ")
  }
  checkmate::reportAssertions(coll)
  rm_df <- colnames(metadata_df)
  rm_df <- rm_df[rm_df != key_col]

  df <- df[,setdiff(colnames(df), rm_df)]
  column_align <- key_col
  names(column_align) <- flowobj_col

    #return(list("flowobj_col" = flowobj_col, "df" = df, "metadata_df" = metadata_df, "key_col" = key_col))
  #df <- dplyr::left_join(df, metadata_df, by = column_align) %>%
  df <- merge(x = df,  y= metadata_df, by.x = flowobj_col, by.y = key_col, all.x = T) %>%
                  tibble::column_to_rownames("tmp_sampleid")

  na_count <- length(apply(is.na(df), 2, which))
  if(na_count > 0){
    print("Warning: there are detected missing values in the metadata. Double check if the supplied dataframe was missing some samples.")
  }
  out_object <- flow_object

  flowWorkspace::pData(out_object$flowSet) <- df
  return(out_object)
}


#' Select files for flow object - Shiny
#'
#' Select files
#'
#' @param flow_object working directory
#'
#' @return A file list
#'

add_metadata_shiny <- function(flow_object){

  ui <- fluidPage(
    headerPanel(
      "File selection"
    ),
    sidebarLayout(
      sidebarPanel(
        tags$h5("Select metadata file"),

        fileInput(inputId = "file", label = "File Select", multiple = F, accept = c("xls","xlsx", "tsv", "csv")),
        # shinyFiles::shinyFilesButton("file", "File select", "Please select a file",
        #                              multiple = FALSE,
        #                              filetype = list(data = c("xls","xlsx", "tsv", "csv")),
        #                              viewtype = "detail"),
        helpText("Note: You can select XLS, XLSX, TSV or CSV files  \n\n"),

        uiOutput("column_sel"),
        uiOutput("column_sel_flow"),
        # tags$h5(htmlOutput("column_valid")),

        actionButton("submit", "Submit")

      ),
      mainPanel(
        tags$h4("Sample Metadata"),
        DT::dataTableOutput("metaDF"),

        tags$hr(),
        tags$h4("Current flow object metadata"),
        DT::dataTableOutput("metaDF_flow")
      )
    )
  )
  server <- function(input, output){
    # volumes <- c(workDir = getwd(), Home = fs::path_home(), shinyFiles::getVolumes()())
    # shinyFiles::shinyFileChoose(input, "file", roots = volumes)

    data <- eventReactive(input$file, {
      #pathfile <- as.character(shinyFiles::parseFilePaths(volumes, input$file)$datapath)
      pathfile <- as.character(input$file$datapath)
      if(grepl("xlsx$|xls$", pathfile, ignore.case = T)){
        df <- readxl::read_excel(pathfile)
      } else if(grepl("tsv$", pathfile, ignore.case = T)){
        df <- readr::read_tsv(pathfile)
      } else {
        df <- readr::read_csv(pathfile)
      }

    })

    metaF <- flowCore::pData(flow_object$flowSet)
    output$metaDF <- DT::renderDataTable({  DT::datatable(data(),
                                                          editable = T,
                                                          options = list(pageLength = 10))})


    output$metaDF_flow <- DT::renderDataTable({  DT::datatable(metaF ,
                                                          editable = F,
                                                          options = list(pageLength = 10))})


    output$column_sel_flow = renderUI({
      if(!is.null(data())){
        selectizeInput("column_sel_flow",
                       label = h5("Select Annotation Column key from FlowSet"),
                       choices = unique(c("name", colnames(as.data.frame(metaF)))),
                       multiple = F,
                       selected = "name",
                       options = list(maxItems = 1))
      }
    })

    output$column_sel = renderUI({
      if(!is.null(data())){
        selectizeInput("column_sel",
                       label = h5("Select Annotation Column"),
                       choices = unique(c("name", colnames(as.data.frame(data())))),
                       multiple = F,
                       selected = "name",
                       options = list(maxItems = 1))
      } else {
        selectizeInput("column_sel",
                       label = h5("Select Annotation Column"),
                       choices = c("name"),
                       multiple = F,
                       selected = "name",
                       options = list(maxItems = 1))
      }
    })

    # output$column_valid <- renderText({
    #   pData_cols <- flowCore::pData(flow_object$flowSet)
    #   if(length(intersect(input$column_sel, colnames(as.data.frame(data())))) == 0){
    #     print( paste("<B>WARNING</B>: Metadata contains no column ", input$column_sel,
    #                  " Select a valid column. \n\n", sep = ""))
    #   } else if(length(intersect(pData_cols$name, as.data.frame(data())[,input$column_sel])) == length(pData_cols$name)) {
    #     print(paste("Matching values found between metadata column ", sQuote(input$column_sel), " and flowSet sample names. Joining...", sep  =""))
    #   } else if(length(intersect(pData_cols$name,
    #                              as.data.frame(data())[,input$column_sel])) == length(pData_cols$name)) {
    #     print(paste("Matching values found between metadata column ", sQuote(input$column_sel), " and flowSet filenames. Joining...", sep = ""))
    #   } else {
    #     print("<B>WARNING</B>: Incomplete or absent matching values from current flowSet metadata. \n\n")
    #   }
    #
    # })


    submitInput <- observeEvent( input$submit,{
      out_df <- as.data.frame(data())
      out <- list("metadata_df" = out_df, "key" = input$column_sel, "fS_key" = input$column_sel_flow )
      shiny::stopApp(out)

    })
  }

  sel <- shiny::runApp(shinyApp(ui = ui, server = server))

}

