#' @import checkmate
#' @import flowCore
#' @import magrittr
#' @import shiny
#' @importFrom rlang .data
#'
NULL



#' Select files for flow object - Shiny
#'
#' Select files
#'
#' @param dir working directory
#'
#' @return A file list
#'

select_files_shiny <- function(dir){

  ui <- fluidPage(
    headerPanel(
      "File selection"
    ),
    sidebarLayout(
      sidebarPanel(
        tags$h5("Select sample files"),


        shinyFiles::shinyFilesButton("file", "File select", "Please select a file",
                                     multiple = TRUE,
                                     filetype = list(data = c("fcs","csv")),
                                     viewtype = "detail"),
        helpText("Note: You can select FCS or CSV files (if exported from FlowJo) \n\n"),

        actionButton("submit", "Submit")

      ),
      mainPanel(
        tags$h4("Selected files"),
        DT::dataTableOutput("fileDF"),

        tags$hr()
      )
    )
  )
  server <- function(input, output){
    volumes <- c(workDir = dir, Home = fs::path_home(), shinyFiles::getVolumes()())
    shinyFiles::shinyFileChoose(input, "file", roots = volumes)

    data <- eventReactive(input$file, {
      pathfile <- as.character(shinyFiles::parseFilePaths(volumes, input$file)$datapath)
      namefile <- as.character(shinyFiles::parseFilePaths(volumes, input$file)$name)
      df <- data.frame("Sample Name" = namefile, "Path" = pathfile)
    })
    output$fileDF <- DT::renderDataTable({  DT::datatable(data(),
                                                              editable = T,
                                                              options = list(pageLength = 30))})
    submitInput <- observeEvent( input$submit,{
      out_df <- as.data.frame(data())
      shiny::stopApp(out_df$Path)

    })
  }

  sel <- shiny::runApp(shinyApp(ui = ui, server = server))

}



#' Create a FlowObject from flow data
#'
#' Randomly selects n events from each file in either a flowSet or a list of FCS files.
#'
#' @param flowSet A FlowSet object. Defaults to NULL. Either a flowSet or files must be supplied using the `files` parameter.
#' @param files A list of flow cytometry sample files. Can be either CSV or FCS files exported by FlowJo. Either a flowSet or files must be supplied using the `flowSet` parameter.
#' @param name identifier for the dataset
#' @param subsample Number of events to subsample per sample. Defaults to 5000.
#' @param min_events Minimum number of events per file. Any sample found with events under this value will be ignored. Defaults to 500.
#' @param max_total_events Total number of events across samples cannot exceed this value: subsampling number will be adjusted accordingly. Defaults to  2e+05
#' @param unequal Logical argument to determine whether or not to allow to take unequal number of events if some samples are limiting. When set to True. Defaults to FALSE.
#' @param interactive Logical argument to determine whether or not to use the interactive mode to select files in a graphical user interface . Defaults to FALSE.
#'
#' @return A Flow Object
#'
#' @examples
#' library(flowCore)
#' data(GvHD)
#' flow_obj <- CreateFlowObject(flowSet = GvHD, name = "GvHD")
#'
#' @export
#'
#'
CreateFlowObject <- function(flowSet = NULL,
                             files = NULL,
                             name = NULL,
                             subsample = 5000,
                             min_events = 500,
                             max_total_events = 200000,
                             unequal = FALSE,
                             interactive = FALSE){

  coll = checkmate::makeAssertCollection()
  checkmate::assertCharacter(name, len = 1, any.missing = F, .var.name = "name", add = coll )
  #assertCharacter(exclude_markers, any.missing = F, .var.name = "exclude_markers", add = coll )

  checkmate::assertNumeric(subsample, len = 1, lower = 100, upper = 1000000, all.missing = F, .var.name = "subsample", add = coll)
  checkmate::assertNumeric(min_events, len = 1, lower = 100, upper = 100000, all.missing = F, .var.name = "min_events", add = coll)
  checkmate::assertCharacter(files, min.len = 1,unique = T, null.ok = T, any.missing = F, .var.name = "files", add = coll)
  checkmate::assertNumeric(max_total_events, len = 1, lower = 100, upper = 20000000, all.missing = F, .var.name = "max_total_events", add = coll)
  checkmate::assertLogical(unequal, any.missing = F, null.ok = F, .var.name = "unequal", add = coll)
  checkmate::assertLogical(interactive, any.missing = F, null.ok = F, .var.name = "interactive", add = coll)
  if(is.null(flowSet) & is.null(files) & interactive == FALSE){
    # print("Insufficient valid files (< 3) in input directory with enough events. Change subsampling parameter or verify files... ABORTING. ")
    coll$push("Error: Must supply at least a value to the 'files'OR 'flowSet' parameters, or use the interactive mode.")
  }
  if(!is.null(flowSet) & !is.null(files)){
    # print("Insufficient valid files (< 3) in input directory with enough events. Change subsampling parameter or verify files... ABORTING. ")
    coll$push("Error: Must supply a value to the 'files' OR 'flowSet' parameters.")
  }


  checkmate::reportAssertions(coll)

  if(!is.null(files) | interactive == TRUE){
    if(interactive == TRUE){
      file_in <- select_files_shiny(getwd())
    } else {
      file_in <- files

    }
    checkmate::assertFile(file_in, .var.name = "files", add = coll)
    checkmate::reportAssertions(coll)
    exten <- tolower(unique(tools::file_ext(file_in)))
    if(length(intersect(exten, c("fcs", "csv"))) != 1){
      coll$push("Error: 'files' can only contain file.paths for either 'csv' OR 'fcs' files.")
    }
    checkmate::reportAssertions(coll)
    if(exten == "csv"){
      flowSet_inp <- csv_to_flowSet(file_in)
    } else {
      flowSet_inp <- flowCore::read.flowSet(files = file_in, transformation = FALSE, truncate_max_range = FALSE)
    }
  }

  if(!is.null(flowSet)){
    if(methods::is(flowSet) != "flowSet"){
      # print("Insufficient valid files (< 3) in input directory with enough events. Change subsampling parameter or verify files... ABORTING. ")
      coll$push(paste("Error: flowSet is of class `", methods::is(flowSet), "`. Supply a valid flowSet generated by `read.flowSet()`", sep = ""))
    }
    checkmate::reportAssertions(coll)

    flowSet_inp <- flowSet
  }

  filenames <- flowCore::sampleNames(flowSet_inp)
  flowCore::pData(flowSet_inp)$filename <- basename(filenames)
  flowCore::pData(flowSet_inp)$name <- paste("Sample", seq(1, length(filenames)), sep = "_")
  flowCore::sampleNames(flowSet_inp) <- paste("Sample", seq(1, length(filenames)), sep = "_")


  minLines <- subsample


  fileLS <- flowWorkspace::sampleNames(flowSet_inp)


  if((length(fileLS) * subsample) > max_total_events) {
    events <- floor((max_total_events/length(fileLS))/50) * 50  ### avoid weird numbers.i.e. subsampling 2941 events per file if the event_number is beyond the max total.
    if(events < minLines){

      print(paste("Subsamping size of ", as.character(subsample),
                  " exceeds the total number of events of ",
                  as.character(max_total_events),
                  ". Rounding down to ", as.character(events), " events per file", sep="" ))
      minLines <- events
    }
  }

  #         # Scanning step: if the fcs files were derived from a pregated population that has a frequency below         # the subsample #, the subsample # needs to be reduced to that number for ALL files. Also used to
  #         # filter out samples with a too small n. of events.

  rmList <- c()

  for (i in fileLS){
    # fcs <- read.FCS(i, truncate_max_range = F)
    fcs <- flowSet_inp[[i]]
    numLines <- dim(fcs@exprs)[1]
    if(numLines < min_events){

      print(paste("WARNING: File ", i, " has ", as.character(numLines),
                  " events, which is lower than the minimum number of events",
                  as.character(min_events), " required. REMOVING FILE FROM ANALYSIS..", sep="" ))
      rmList <- c(rmList,i)

    } else if(unequal == FALSE  & numLines < minLines){
      print(paste("WARNING: File ", i, " has ", as.character(numLines), " events,
                  lower than specified current subsampling number of ", as.character(minLines),
                  ". Changing subsampling parameter...", sep=""))
      minLines <- numLines
    }
  }

  # Filter out files with  insufficient events and check if enough are left
  fileLS_filt <- fileLS[!fileLS %in% rmList]
  if(length(fileLS_filt) < 1){
    coll$push("Insufficient valid files (< 1) in input directory with enough events. Change subsampling parameter ... ABORTING. ")
  }
  checkmate::reportAssertions(coll)


  fS_filt <- Downsampling_FlowSet(flowSet_inp[fileLS_filt], n = minLines )
  for(i in names(fS_filt@frames)){
    row.names(flowCore::exprs(fS_filt[[i]])) <- paste(i, seq(1, to = dim(flowCore::exprs(fS_filt[[i]]))[1]), sep ="__")

  }

  if(suppressWarnings(class(flowCore::markernames(fS_filt))) == "list"){
    coll$push("Marker names are not consistent across samples within flowSet. Aborting...")
    checkmate::reportAssertions(coll)
  }

  #exclude_marker_check <- exclude_markers[!exclude_markers %in% c(unname(markernames(fS_filt)), unname(pData(parameters(flowSet[[1]]))$name))]
  #print(exclude_marker_check)
  # if(length(exclude_marker_check)> 0){
  #   coll$push(paste("Excluded markers ", paste(sQuote(exclude_marker_check), collapse = ", "), " are not in the subsample files.", sep = ""))
  # }
  # reportAssertions(coll)

  param_LS <- list("exclude_markers" = list(),
                   "include_markers" = flowCore::markernames(flowSet_inp),
                   "cofactor_find" = FALSE)

  output <- list("flowSet"= fS_filt,
                 "dataset" = name,
                 "parameters" = param_LS)

  class(output) <- "flow_object"

  return(output)

}


#' Get fluorescence data
#'
#' Extract the expression matrix from a flow_object, which can be complemented with metadata.
#'
#' @param flow_object A flow_object
#' @param add_sample_id Logical argument to determine whether or not to add a SampleID column for each cell
#' @param annotations Logical argument to determine whether or not to add additional annotations to the data, including clusters, cluster annotations and sample metadata, and dimensions reductions.
#'
#' @return A data.frame
#'
#' @examples
#' library(flowCore)
#' data("GvHD")
#' flow_obj <- CreateFlowObject(flowSet = GvHD, name = "GvHD")
#' exprs <- get_flowSet_matrix(flow_obj, add_sample_id = TRUE)
#' @export
#'
#'

get_flowSet_matrix <- function(flow_object,
                               add_sample_id = FALSE,
                               annotations = FALSE){

  coll <- checkmate::makeAssertCollection()
  if(methods::is(flow_object) != "flow_object") { coll$push("flow_object not recognized. Please supply a valid flow panel. ")}
  checkmate::assertLogical(add_sample_id, len = 1, .var.name = "add_sample_id")
  checkmate::assertLogical(annotations, len = 1, .var.name = "annotations")
  checkmate::reportAssertions(coll)
  flowDF <- do.call("rbind", lapply(names(flow_object$flowSet@frames), function(y){
    df <- as.data.frame(flowCore::exprs(flow_object$flowSet[[y]]))
    #if(add_sample_id == TRUE){
    df$SampleID <- rep(y, dim(df)[1])
    #df$cellID <- paste(df$SampleID, row.names(df), sep ="__")
    #}
    return(df)
  })) %>%
    tibble::rownames_to_column("cellID")

  #extra_cols <- c("UMAP1", "UMAP2", "louvain")

  channels <- intersect(colnames(flowDF), names(flowCore::markernames(flow_object$flowSet)))

  markers <- unname(flowCore::markernames(flow_object$flowSet))

  # if(annotations){
  #   add_sample_id <- TRUE
  # }

  #if(add_sample_id){
  channels <- c(channels, "SampleID", "cellID")
  markers <- c(markers, "SampleID", "cellID")
  #}
  flowDF <- flowDF[,channels]
  colnames(flowDF) <- markers


  if(!is.null(flow_object$dims)){
    flowDF <- flowDF %>%
      dplyr::left_join(flow_object$dims %>%
                         tibble::rownames_to_column("cellID"), by = "cellID")
  }
  if(!is.null(flow_object$louvain)){
    flowDF <- flowDF %>%
      dplyr::left_join(flow_object$louvain %>%
                         tibble::rownames_to_column("cellID"), by = "cellID")
  }

  if(annotations == T){
    flowDF <- flowDF %>%
      dplyr::left_join(as.data.frame(flowWorkspace::pData(flow_object$flowSet)) %>%
                         tibble::rownames_to_column("SampleID"),
                       by = "SampleID") %>%
      dplyr::select(-.data$SampleID)
  }
  if(add_sample_id == F){
    flowDF <- flowDF %>% dplyr::select(-.data$SampleID, -.data$cellID)
  }
  return(flowDF)
}


#' Get cluster frequencies
#'
#' Extract cluster frequencies from a flow object. If the object contains cluster annotations, the function will calculate frequencies for the indicated annotation.
#'
#' @param flow_object A flow_object
#' @param louvain_select A filter for a subset of louvain clusters for which you want frequencies. NULL by default.
#' @param annotation Add additional annotations to the data, including clusters, cluster annotations and sample metadata, and dimensions reductions. Set to 'louvain' by default.
#' @param format The format of the output. Can be set to 'long' or 'wide'. Annotations can only be added in a long format. Defaults to 'long'.
#' @param output Retrieve relative or absolute frequencies. If the downsampling was NOT equalized at flow object creation, relative frequencies are recommended.
#'
#' @return A data.frame
#'
#' @export
#'
#'

get_cluster_frequencies <- function(flow_object,
                                    louvain_select = NULL,
                                    annotation = "louvain",
                                    format = "long",
                                    output = "relative"){
  coll <- checkmate::makeAssertCollection()
  if (methods::is(flow_object) != "flow_object") {
    coll$push("flow_object not recognized (nor container or panel). Please supply a valid flow object.")
  }

  checkmate::assertCharacter(as.character(louvain_select), null.ok = T ,.var.name = "louvain_select", add = coll)
  checkmate::assertCharacter(format, null.ok = F, any.missing = F, len = 1 ,.var.name = "format", add = coll)
  checkmate::assertChoice(format, choices = c("long", "wide"), .var.name = "format", add = coll)
  checkmate::assertCharacter(output, null.ok = F, any.missing = F, len = 1 ,.var.name = "output", add = coll)
  checkmate::assertChoice(output, choices = c("absolute", "relative"), .var.name = "output", add = coll)
  checkmate::assertCharacter(annotation, null.ok = F, any.missing = F, len = 1 ,.var.name = "annotation", add = coll)

  checkmate::reportAssertions(coll)

  flowDF_annot <- get_flowSet_matrix(flow_object, add_sample_id = T)

  if(!"louvain" %in% colnames(flowDF_annot)){
    coll$push(paste("Louvain clusters not found in flow object. First run ",sQuote("Cluster_flow"), ". ", sep = ""))
  }

  flowDF_annot <- flowDF_annot[,c("louvain", "SampleID", "cellID")] #%>%
  # dplyr::group_by(SampleID, louvain) %>%
  # dplyr::summarise(n_total = n_distinct(cellID)) %>%
  # dplyr::ungroup() %>%
  # tidyr::spread(SampleID, n_total, fill = 0) %>%
  # tidyr::gather(SampleID, n_total, -louvain)

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
    annotation_valid <- colnames(clust_annot)
    if(!annotation %in% annotation_valid){
      coll$push(paste("Chosen annotation ", sQuote(annotation), " is not in current cluster annotations. ",
                      "Valid values are ", paste(sQuote(annotation_valid), collapse = ";"), sep =""))
    }
    checkmate::reportAssertions(coll)

  }

  group_index <- which(colnames(df) == annotation)
  sample_index <- which(colnames(df) == "SampleID")

  sample_counts <- df %>%
    dplyr::group_by(.data$SampleID) %>%
    dplyr::summarise(n_total = dplyr::n()) %>%
    dplyr::ungroup()
  df <- df %>%
    dplyr::select(dplyr::all_of(c(annotation, "SampleID","cellID"))) %>%
    dplyr::group_by_at(vars(dplyr::all_of(annotation),"SampleID" )) %>%
    dplyr::summarise(absolute_counts = dplyr::n()) %>%
    dplyr::ungroup() %>%
    tidyr::spread(.data$SampleID, .data$absolute_counts, fill = 0) %>%
    tidyr::gather(key = "SampleID", value = "absolute_counts", -1) %>%
    dplyr::left_join(sample_counts, by = "SampleID") %>%
    dplyr::mutate(relative_frequency = .data$absolute_counts / .data$n_total) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.data$n_total) %>%
    base::unique() #%>%

  df <- df %>% dplyr::left_join(sample_metadata[,c(setdiff(colnames(sample_metadata), colnames(df)), "SampleID" )], by = "SampleID")
  #return(df)
  if("louvain" %in% colnames(df)){
    if(!is.null(louvain_select)){
      invalid_filter <- setdiff(as.character(louvain_select), as.character(unique(df$louvain)))
      if(length(invalid_filter) > 0){
        coll$push(paste("Invalid louvain cluster filter values :" ,
                        paste(invalid_filter, collapse = ";"),
                        ". Valid values range from ", sQuote(min(df$louvain)), " to ",
                        sQuote(max(df$louvain)), ". ", sep = ""))
      }
      df <- df[df$louvain %in% as.numeric(louvain_select),]

    }
    df$louvain <- factor(as.character(df$louvain))

  }

  group_reorder <- c("SampleID", annotation, "absolute_counts", "relative_frequency")

  reorder_ls <-  c(which(colnames(df) %in% group_reorder), which(!colnames(df) %in% group_reorder))

  #df <- df %>% dplyr::select(.data$reorder_ls)
  df <- df[,reorder_ls]


  return(df)

}

#' Transform fluorescence values
#'
#' Applies a transformation to fluorescence values embeeded in a Flow Object. Possible choices are hyperbolic sine (asinh)
#' using cofactors inferred by `find_cofactors()`, biexponential (none), logarithmic (in base 10), or CYTOF asinh (cofactor of 5).
#'
#' @param flow_object A flow_object
#' @param transform Transformation method. Possible choices are 'asinh', 'none', 'log' or 'cytofAsinh'.
#'
#' @return A data.frame
#'
#' @export
#'
#'
transform_flow_data <- function(flow_object,
                                transform ){
  coll <- checkmate::makeAssertCollection()
  if(methods::is(flow_object) != "flow_object") { coll$push("flow_object not recognized. Please supply a valid flow panel. ")}
  checkmate::assertCharacter(transform, len = 1, any.missing = F, null.ok = F, .var.name = "transform", add = coll)
  valid_trans <- c("asinh", "none", "log", "cytofAsinh")
  if(!transform %in% valid_trans){
    coll$push(paste("Transformation ", sQuote(transform), " is not a valid value. Valid options are ", paste(valid_trans, collapse = "; "),
                    ".", sep = ""))
  }
  #assertChoice(transform, choices = c("asinh", "none", "logicle", "cytofAsinh", "log"), .var.name = "transform", add = coll )
  checkmate::reportAssertions(coll)

  flowDF_annot <- get_flowSet_matrix(flow_object, add_sample_id = T)
  include <- flow_object$parameters$include_markers
  flowDF <- flowDF_annot[,include]

  if(transform == "asinh"){
    if(flow_object$parameters$cofactor_find == FALSE){
      coll$push("No cofactors found in FlowObject for hyperbolic sine transformation. First run `find_cofactors()`")
      checkmate::reportAssertions(coll)
    } else{
      cofactors <- flow_object$parameters$cofactors
    }
    flowDF <- do.call("cbind", lapply(flow_object$parameters$include_markers,function(x){
      df <- data.frame(asinh(flowDF[,x]/cofactors[[x]]))
      colnames(df) <- x
      return(df)
    }))

  } else if(transform == "cytofAsinh"){
    cofactors <- rep(5, length(include))
    names(cofactors) <- include
    flowDF <- do.call("cbind", lapply(names(cofactors),function(x){
      df <- data.frame(asinh(flowDF[,x]/cofactors[[x]]))
      colnames(df) <- x
      return(df)
    }))
  } else if(transform == "log"){
    ### function from CytoTree
    ### https://github.com/JhuangLab/CytoTree/blob/master/R/preprocessing.R
    cs <- apply(abs(flowDF[ ,include]), 2, sum)
    #transMarker_id <- transMarker_id[cs > 0]
    #cs <- cs[cs > 0]
    norm_factors <- (10**ceiling(log10(stats::median(cs))))/cs
    flowDF <- round(log10(base::sweep(abs(flowDF[ ,include]), 2, norm_factors, "*")+1), digits=3)
  }

  row.names(flowDF) <- flowDF_annot$cellID
  return(flowDF)
}


#' Exclude markers
#'
#' Excludes markers in order to avoid using them for dimension reduction, clustering and other techniques.
#' To be used if cells were pregated on, or if some markers are unreliable.
#'
#' @param flow_object A Flow Object
#' @param exclude_markers A vector of markernames to exclude.
#'
#' @return A Flow Object
#'
#' @export
#'

exclude_parameters <- function(flow_object,
                               exclude_markers) {
  coll <- checkmate::makeAssertCollection()
  if(methods::is(flow_object) != "flow_object") { coll$push("flow_object not recognized. Please supply a valid flow panel. ")}
  checkmate::reportAssertions(coll)

  invalid_params <- exclude_markers[!exclude_markers %in% flowCore::markernames(flow_object$flowSet)]

  if(length(invalid_params)>0){
    coll$push(paste("Parameter names ", paste(sQuote(invalid_params), collapse = "; "), " are not valid parameters in the flow Object.", sep = ""))
  }

  checkmate::reportAssertions(coll)
  flow_object$parameters$exclude_markers <- unlist(unique(c(flow_object$parameters$exclude_markers,exclude_markers)))
  flow_object$parameters$include_markers <- setdiff(flowCore::markernames(flow_object$flowSet), exclude_markers)


  return(flow_object)

}


