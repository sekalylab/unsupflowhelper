#' @include generics.R
#' @import magrittr
#' @importFrom rlang .data
#'
NULL

#' Regression analysis of cluster frequencies
#'
#' Statistical comparison of cluster frequencies versus metadata continuous variables. Metadata needs to be added using the `add_sample_metadata` function.
#'
#' @param flow_object A Flow Object
#' @param outcome Numerical metadata variable or outcome
#' @param group.by Which cluster annotation to compute with: defaults to 'louvain'
#' @param method Correlation test to use: either 'spearman' or 'pearson'. Defaults to `spearman`.
#'
#' @return A data.frame
#'
#' @export
#'

correlate_clusters <- function(flow_object,
                               outcome,
                               method = "spearman",
                               group.by = "louvain"
                               ){

  coll <- checkmate::makeAssertCollection()
  if (methods::is(flow_object) != "flow_object") {
    coll$push("flow_object not recognized (nor container or panel). Please supply a valid flow object.")
  }

  checkmate::assertCharacter(method, len = 1, any.missing = F, null.ok = F,.var.name = "method", add = coll )
  checkmate::assertCharacter(group.by, len = 1, any.missing = F, null.ok = F, .var.name = "group.by", add = coll)
  checkmate::assertChoice(method, choices = c("spearman", "pearson"), .var.name = "method", add = coll)
  checkmate::assertCharacter(outcome, len = 1, any.missing = F, null.ok = F,.var.name = "outcome", add = coll )
  checkmate::reportAssertions(coll)

  df <- get_cluster_frequencies(flow_object, annotation = group.by)

  numeric_columns <- colnames(df)[unlist(lapply(df, is.numeric), use.names = FALSE)]
  character_column <- colnames(df)[!unlist(lapply(df, is.numeric), use.names = FALSE)]
  numeric_columns <- setdiff(numeric_columns, c("absolute_counts", "relative_frequency"))

  if(!outcome %in% numeric_columns ){
    coll$push(paste("Variable ", outcome, " is not a valid numerical metadata variable", sep =""))
  } else if(!group.by %in% character_column){
    coll$push(paste("Variable ", group.by, " is not a valid discrete grouping metadata variable", sep =""))
  }
  checkmate::reportAssertions(coll)


  df <- df %>%
    tidyr::gather(key = "outcome_var", value = "outcome_value", dplyr::all_of(outcome)) %>%
    tidyr::gather(key = "grouping", value = "grouping_value", dplyr::all_of(group.by)) %>%
    dplyr::group_by(.data$grouping_value, .data$SampleID) %>%
    dplyr::mutate(sum_freq = sum(.data$relative_frequency)) %>%
    dplyr::select(-.data$relative_frequency, -.data$absolute_counts) %>%
    dplyr::ungroup()

  df <- df %>%
    dplyr::group_by(.data$grouping_value) %>%
    rstatix::cor_test(.data$outcome_value, .data$sum_freq, method = method ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(padj = stats::p.adjust(.data$p, method = "BH")) %>%
    dplyr::rename(p.value = "p") %>%
    dplyr::select(-.data$var2) %>%
    dplyr::mutate(var1 = outcome) %>%
    dplyr::rename(outcome = "var1") #%>%
  colnames(df)[1] <- group.by

  if(method == "spearman"){
    df <- df %>% dplyr::rename(rho = "cor")

  } else {
    df <- df %>% dplyr::rename(p = "cor")
  }



  return(df)
}


#' Discrete analysis of cluster frequencies
#'
#' Statistical comparison of cluster frequencies versus metadata discrete (i.e. categorical) variables. Metadata needs to be added using the `add_sample_metadata` function.
#'
#' @param flow_object A Flow Object
#' @param outcome Discrete or categorical metadata variable or outcome. Can be of 'character' or 'factor' class.
#' @param group.by Which cluster annotation to compute with: defaults to 'louvain'
#' @param method Correlation test to use: either 'wilcox' or 't.test'. Defaults to `wilcox`. Columns with more than 2 levels of categories will use the 'pairwise' versions fo the statistical tests.
#'
#' @return A data.frame
#'
#' @export
#'


clusters_compare_groups <- function(flow_object,
                                    outcome,
                                    method = "non-parametric",
                                    group.by = "louvain"
                                    ){
  coll <- checkmate::makeAssertCollection()
  if (methods::is(flow_object) != "flow_object") {
    coll$push("flow_object not recognized (nor container or panel). Please supply a valid flow object.")
  }

  checkmate::assertCharacter(method, len = 1, any.missing = F, null.ok = F,.var.name = "method", add = coll )
  checkmate::assertCharacter(group.by, len = 1, any.missing = F, null.ok = F, .var.name = "group.by", add = coll)
  checkmate::assertChoice(method, choices = c("parametric", "non-parametric"), .var.name = "method", add = coll)
  checkmate::assertCharacter(outcome, len = 1, any.missing = F, null.ok = F,.var.name = "outcome", add = coll )
  checkmate::reportAssertions(coll)

  df <- get_cluster_frequencies(flow_object, annotation = group.by)

  numeric_columns <- colnames(df)[unlist(lapply(df, is.numeric), use.names = FALSE)]
  character_columns <- colnames(df)[!unlist(lapply(df, is.numeric), use.names = FALSE)]
  numeric_columns <- setdiff(numeric_columns, c("absolute_counts", "relative_frequency"))

  if(outcome %in%  numeric_columns ){
    coll$push(paste("Variable ", outcome, "is not a valid discrete metadata variable", sep =""))
  } else if(!group.by %in% character_columns){
    coll$push(paste("Variable ", group.by, "is not a valid discrete grouping metadata variable", sep =""))
  }
  checkmate::reportAssertions(coll)

  if(is.factor(df[[outcome]])){
    fact_levs <- levels(df[[outcome]])
  } else {
    fact_levs <- levels(factor(df[[outcome]]))
  }

  group_n <- as.data.frame(table(flowWorkspace::pData(flow_object$flowSet)[outcome]))

  if(min(group_n$Freq) < 3){
    small_groups <- unique(as.character(group_n[group_n$Freq < 3,][[1]]))
    coll$push(paste("Error: Group(s):", paste(sQuote(small_groups), collapse = "; "), "for outcome ",
                    sQuote(outcome), " have fewer than 3 samples."))
  }
  checkmate::reportAssertions(coll)

  df <- df %>%
    tidyr::gather(key = "outcome_var", value = "outcome_value", dplyr::all_of(outcome)) %>%
    dplyr::mutate(outcome_value = factor(.data$outcome_value, levels = fact_levs)) %>%
    tidyr::gather(key = "grouping", value = "grouping_value", dplyr::all_of(group.by)) %>%
    dplyr::group_by(.data$grouping_value, .data$SampleID) %>%
    dplyr::mutate(sum_freq = sum(.data$relative_frequency)) %>%
    dplyr::select(-.data$relative_frequency, -.data$absolute_counts) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$grouping_value) %>%
    dplyr::mutate(n_samples= dplyr::n_distinct(.data$SampleID))


  ngroups <- length(unique(df[["outcome_value"]]))


  if(ngroups == 1){
    coll$push("Outcome variable only contains one level and cannot be used for contrasts")
    checkmate::reportAssertions(coll)
  } else if(ngroups >2){

    if(method == "parametric"){
      df <- df %>%
        rstatix::t_test(stats::as.formula("sum_freq ~ outcome_value"), paired = F, data = .data) %>%
        #tukey_hsd(aov(sum_freq ~ outcome_value, data = .)) %>%
        dplyr::mutate(Test = "One way Anova with Tukey HSD post-hoc test")
    } else{
      df <- df %>%
        rstatix::wilcox_test( stats::as.formula("sum_freq ~ outcome_value"), data = .data) %>%
        dplyr::mutate(Test = "Non-parametric pairwise Wilcoxon-tests")

    }
  } else if(ngroups == 2){
    if(method == "parametric"){
      df <- df %>%
        rstatix::t_test(stats::as.formula("sum_freq ~ outcome_value")) %>%
        dplyr::mutate(Test = "Parametric 2 tailed T-test")
    } else{
      df <- df %>%
        rstatix::wilcox_test(stats::as.formula("sum_freq ~ outcome_value")) %>%
        dplyr::mutate(Test = "Non-parametric Wilcoxon Test")
    }

  }

  df <- df %>%
    dplyr::ungroup() %>%
    dplyr::rename(p.value = "p") %>%
    dplyr::select(-.data$`.y.`) %>%
    #mutate(var1 = outcome) %>%
    dplyr::mutate(outcome = outcome) #%>%
  colnames(df)[1] <- group.by


  return(df)
}
