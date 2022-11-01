#' @include generics.R
#' @include utilities.R
#' @import flowCore
#' @import magrittr
#' @import shiny
#' @importFrom rlang .data
#'
NULL


#' Export a Flow Object to a FCS file
#'
#' Generates a single FCS file, that can be imported into FlowJo, with embedded cluster labels, UMAP coordinates, and sample metadata.
#' Columns with character values, like cluster annotations or sample annotations, will be encoded as numbers, and a key to those values will be generated as a seaparate CSV file.
#'
#' @param flow_object A Flow Object
#' @param filename A filename for the output FCS file. Must contain the `.fcs` extension.
#'
#' @return A FCS file
#'
#' @export
#'

export_to_fcs <- function(flow_object,
                          filename){

  coll <- checkmate::makeAssertCollection()

  if(methods::is(flow_object) != "flow_object"){
    coll$push("flow_object not recognized. Please supply a valid flow object.")
  }

  checkmate::assertCharacter(filename, add = coll, len = 1, any.missing = F, null.ok = F, .var.name = "filename")
  if(grepl(".fcs$", filename, ignore.case = T) == FALSE){
    coll$push("Invalid Filename: must have `.fcs` entension")
  }
  checkmate::reportAssertions(coll)

  fS <- flow_object$flowSet
  fS_names <- colnames(fS)

  flowDF_annot <- get_flowSet_matrix(flow_object, add_sample_id = T, annotations = F)
  if("louvain" %in% colnames(flowDF_annot)){
    flowDF_annot$louvain <- as.numeric(flowDF_annot$louvain)
  }
  flowDF_annot <- flowDF_annot[,intersect(c("UMAP1", "UMAP2", "louvain", "SampleID"), c("cellID",colnames(flowDF_annot)))]


  sample_metadata <- pData(flow_object$flowSet) %>% tibble::rownames_to_column("SampleID")

  if(length(setdiff(colnames(sample_metadata), c("SampleID", "name"))) > 0){
    numeric_columns <- colnames(sample_metadata)[unlist(lapply(sample_metadata, is.numeric), use.names = FALSE)]
    character_columns <- colnames(sample_metadata)[unlist(lapply(sample_metadata, is.character), use.names = FALSE)]
    factor_columns <- colnames(sample_metadata)[unlist(lapply(sample_metadata, is.factor), use.names = FALSE)]
    character_columns <- character_columns[!character_columns %in% c("SampleID")]

    factorLS <- list()
    for(i in character_columns){
      sample_metadata[,i] <- factor(sample_metadata[,i])
      factorLS[[i]] <- data.frame("value" = levels(sample_metadata[,i]),
                                  "key" = seq(length(levels(sample_metadata[,i])))) %>%
        dplyr::mutate(metadata_column = i)
      sample_metadata[,i] <- as.numeric(sample_metadata[,i])
    }
    for(i in factor_columns){
      #sample_metadata[,i] <- factor(sample_metadata[,i])
      factorLS[[i]] <- data.frame("value" = levels(sample_metadata[,i]),
                                  "key" = seq(length(levels(sample_metadata[,i])))) %>%
        dplyr::mutate(metadata_column = i)
      sample_metadata[,i] <- as.numeric(sample_metadata[,i])
    }

    key_df <- dplyr::bind_rows(factorLS)

    keyfile <- paste(filename,".key.csv", sep ="")
    utils::write.csv(key_df, file = keyfile, row.names = F)

    flowDF_annot <- flowDF_annot %>%
                    dplyr::left_join(sample_metadata, by = "SampleID")
    factor_mat <- flowDF_annot[,colnames(sample_metadata)] %>%
                  dplyr::select(-.data$SampleID)
    #out_cols <- intersect(c("UMAP1", "UMAP2", "louvain", colnames(sample_metadata)), )
  } else {
    factor_mat <- flowDF_annot %>% dplyr::select(-.data$SampleID)
  }

  fS <- flow_object$flowSet
  fS_names <- colnames(fS)

  fS <- flowCore::flowSet_to_list(fS)

  for(i in names(fS)){
    ff <- fS[[i]]
    rown <- row.names(flowDF_annot[flowDF_annot$SampleID == i,])
    new_frame <- flowCore::fr_append_cols(ff, as.matrix(factor_mat[rown,]))
    fS[[i]] <- new_frame
  }

  fS_new <- methods::as(fS, "flowSet")
  flowWorkspace::pData(fS_new) <- flowWorkspace::pData(flow_object$flowSet)

  ff <- MetaCyto::set2Frame(fS_new)

  limit_reset_cols <- c(c("UMAP1", "UMAP2", "louvain"), colnames(factor_mat))

  param <- as.data.frame(ff@parameters@data) %>%
          tibble::rownames_to_column("rown") %>%
          dplyr::mutate(minRange = ifelse(.data$name %in% limit_reset_cols, .data$minRange -1, .data$minRange),
                        maxRange = ifelse(.data$name %in% limit_reset_cols, .data$maxRange +1, .data$maxRange),
                        range = ifelse(.data$name %in% limit_reset_cols, .data$range + 2, .data$range)) %>%
          tibble::column_to_rownames("rown")
  ff@parameters@data <- param

  flowCore::write.FCS(ff, filename)



  #}


}


#' Export a Flow Object to a spreadsheet
#'
#' Generates a single FCS file, that can be imported into FlowJo, with embedded cluster labels, UMAP coordinates, and sample metadata.
#' Columns with character values, like cluster annotations or sample annotations, will be encoded as numbers, and a key to those values will be generated as a seaparate CSV file.
#'
#' @param flow_object A Flow Object
#' @param filename A filename for the output XLSX file. Must contain the `.xls` or `xlsx` extension.
#'
#' @return A spreadsheet file
#'
#' @export
#'

export_to_spreadsheet <- function(flow_object,
                                  filename){

  coll <- checkmate::makeAssertCollection()

  if(methods::is(flow_object) != "flow_object"){
    coll$push("flow_object not recognized. Please supply a valid flow object.")
  }
  checkmate::assertCharacter(filename, add = coll, len = 1, any.missing = F, null.ok = F, .var.name = "filename")
  if(grepl(".xls(x|)$", filename, ignore.case = T) == FALSE){
    coll$push("Invalid Filename: must have `.xls` or `.xlsx` entension")
  }
  checkmate::reportAssertions(coll)


  flowDF_annot <- get_flowSet_matrix(flow_object, add_sample_id = T)
  marker_cols <- which(colnames(flowDF_annot) %in% setdiff(flowCore::markernames(flow_object$flowSet),
                                                           c("UMAP1", "UMAP2", "louvain", "partition")))
  sample_metadata <- flowWorkspace::pData(flow_object$flowSet) %>% tibble::rownames_to_column("SampleID")
  df <- flowDF_annot %>%
    dplyr::left_join(sample_metadata, by = "SampleID") %>%
    dplyr::ungroup()



  if("louvain" %in% colnames(df)){
    df$louvain <- as.character(df$louvain)
    clust_annot <- flow_object$parameters$louvain_annotations
    if(dim(clust_annot)[2] > 1){
      df <- df %>%
        dplyr::mutate(louvain = as.character(.data$louvain)) %>%
        dplyr::left_join(clust_annot, by = "louvain")

    }
  }
  exportLS <- list("flowDF" = df)
  exportLS$Cluster_counts <- get_cluster_frequencies(flow_object = flow_object)
  exportLS$cluster_MFI <- df %>%
    tidyr::gather(key = "marker", value = "expression", dplyr::all_of(marker_cols)) %>%
    dplyr::group_by(.data$louvain, .data$marker) %>%
    dplyr::summarise(median = stats::median(.data$expression)) %>%
    dplyr::ungroup() %>%
    tidyr::spread(.data$marker, .data$median)

  writexl::write_xlsx(exportLS, path = filename)
}
