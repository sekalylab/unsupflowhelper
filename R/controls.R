#' @include generics.R
#' @import flowCore
#' @import magrittr
#' @import shiny
#' @importFrom rlang .data
#'
NULL

#' Add Staining controls
#'
#' Select control FCS files as references for plotting.
#'
#' @param flow_object A Flow Object
#' @param interactive Activates the interactive mode: opens a graphical interface to select samples instead of providing a control_list. Defaults to TRUE.
#' @param control_list A named list where names correspond to markers found in the Flow Object, and values are paths to FCS files. At least one value must also be for an Unstained Universal Control. see Details.
#' @param downsample Number of events to select from each FCS file.
#'
#' @return A FlowSet
#'
#' @export
#'

add_staining_controls <- function(flow_object,
                                  interactive = T,
                                  control_list = NULL,
                                  downsample = 2000
                                  ){
  coll <- checkmate::makeAssertCollection()
  if(methods::is(flow_object) != "flow_object"){
    coll$push("flow_object not recognized. Please supply a valid flow object.")
  }
  checkmate::assertList(control_list, any.missing = F, null.ok = T, .var.name = "control_list", add = coll)
  checkmate::assertLogical(interactive, any.missing = F, .var.name = "interactive", add = coll)
  checkmate::reportAssertions(coll)
  valid_names <- flowCore::markernames(flow_object$flowSet)
  names(valid_names) <- NULL
  valid_names <- valid_names[!valid_names %in% c("UMAP1", "UMAP2", "louvain")]
  #valid_names <- setdiff(valid_names, flow_object$parameters$exclude_markers)
  if(interactive == TRUE){
    control_LS <- add_staining_controls_shiny(valid_names)

  } else{
    invalid_names <- setdiff(names(control_list), c("Unstained", valid_names))

    if(length(invalid_names) > 0){
      coll$push(paste("Invalid marker names for values: {", paste(invalid_names, collapse = ";"), "}.   ",
                      "Valid values can be {", paste(c("Unstained", valid_names), collapse = "; "), "}."))
    }
    checkmate::reportAssertions(coll)
    if(!"Unstained" %in% names(control_list)){
      coll$push(paste("Reference controls must contain an `Unstained` control. "))
    }
    checkmate::reportAssertions(coll)
    control_LS <- control_list
  }



  FS <- suppressWarnings(read.flowSet(files= unlist(control_LS), transformation = FALSE, truncate_max_range = FALSE))
  FS <- Downsampling_FlowSet(FS, n = downsample)
  flowWorkspace::sampleNames(FS) = names(control_LS)
  flowWorkspace::pData(FS)$name <- names(control_LS)
  out_object <- flow_object
  out_object$staining_controls <- FS
  return(out_object)
}

#' Get reference matrix
#'
#' Extract fluorescence values from the reference control FlowSet.
#'
#' @param flow_object A Flow Object
#' @param add_sample_id Add a SampleID column to the matrix. Defaults to TRUE.
#'
#' @return A data.frame or matrix.
#'
#' @export
#'

get_reference_matrix <- function(flow_object,
                                 add_sample_id = TRUE){
  coll <- checkmate::makeAssertCollection()
  if(methods::is(flow_object) != "flow_object") { coll$push("flow_object not recognized. Please supply a valid flow panel. ")}
  checkmate::assertLogical(add_sample_id, len = 1, .var.name = "add_sample_id")
  checkmate::reportAssertions(coll)
  flowDF <- do.call("rbind", lapply(names(flow_object$staining_controls@frames), function(y){
    df <- as.data.frame(flowCore::exprs(flow_object$staining_controls[[y]]))
    if(add_sample_id == TRUE){
      df$SampleID <- rep(y, dim(df)[1])
    }
    return(df)
  }))


  channels <- intersect(colnames(flowDF), c(names(flowCore::markernames(flow_object$staining_controls))))
  markers <- unname(flowCore::markernames(flow_object$staining_controls))

  if(add_sample_id){
    channels <- c(channels, "SampleID")
    markers <- c(markers, "SampleID")
  }
  flowDF <- flowDF[,channels]
  colnames(flowDF) <- markers
  return(flowDF)
}

#' Shiny interface to add reference controls
#'
#' Boot the GUI file selector
#'
#' @param markers A Flow Object
#'
#' @return A named list of reference control files.
#'

add_staining_controls_shiny <- function(markers ){
  control_df <- data.frame("Marker" = c("Unstained",sort(markers))) %>%
    dplyr::mutate(name = "",
                  filepath = "")

  ui <- fluidPage(
    headerPanel(
      "Control file selection"
    ),
    sidebarLayout(
      sidebarPanel(
        tags$h5("Select a control and its associated fcs file"),
        selectInput("marker", "Control:",
                    choices = control_df$Marker),

        shinyFiles::shinyFilesButton("file", "File select", "Please select a file",
                         multiple = FALSE,
                         filetype = "fcs",
                         viewtype = "detail"),
        helpText("Note: this functionality requires at least a universal control (i.e. Unstained).",
                 "If additional controls are provided for other markers, the MFI of that marker from that control will be used.",
                 "Controls are not necessary for every marker. \n\n"),

        actionButton("submit", "Submit")

      ),
      mainPanel(
        tags$h4("Selected reference controls"),
        tags$h5(htmlOutput("unstain_valid")),
        DT::dataTableOutput("control_df"),

        tags$hr()
      )
    )
  )


  server <- function(input, output){
    volumes <- c(workDir = getwd(), Home = fs::path_home(), shinyFiles::getVolumes()())
    shinyFiles::shinyFileChoose(input, "file", roots = volumes)

    rv <- reactiveVal(control_df)

    observeEvent(input$file, {
      pathfile <- as.character(shinyFiles::parseFilePaths(volumes, input$file)$datapath)
      namefile <- as.character(shinyFiles::parseFilePaths(volumes, input$file)$name)
      cat("Adding", pathfile, "for", input$marker, "\n")
      new_df <- rv()
      new_df[which(new_df$Marker == input$marker),][2] <- namefile
      new_df[which(new_df$Marker == input$marker),][3] <- pathfile
      rv(new_df)
      updateNumericInput(inputId = "file", value = 1)
    })

    output$control_df <- DT::renderDataTable({  DT::datatable(rv(),
                                                              editable = T,
                                                              options = list(pageLength = 30))})

    output$unstain_valid <- renderText({
      df_unst <- as.data.frame(rv())
      df_unst <- df_unst[df_unst$Marker == "Unstained",]
      if(df_unst$name[1] == ""){
        print("<B>WARNING</B>: `Unstained` control has no selected file. A standard universal control is <B>necessary</B> for this functionality.\n\n")
      } else {
        print("Valid `Unstained` control detected. Click <B>Submit</B> when done selecting control files.\n\n")
      }
    })
    submitInput <- observeEvent( input$submit,{
      out_df <- as.data.frame(rv())
      out_df <- out_df[out_df$name != "",]
      df_unst <- out_df[out_df$Marker == "Unstained",]
      if(df_unst$name != ""){
        outLS <- as.list(out_df$filepath)
        #print(outLS)
        names(outLS) <- out_df$Marker
        #print(outLS)
        #vals <- shiny::reactiveValues(outLS = outLS)
        if(df_unst$name[1] != ""){
          shiny::stopApp(outLS)
        }
      }
    })
  }

  sel <- shiny::runApp(shinyApp(ui = ui, server = server))

}

