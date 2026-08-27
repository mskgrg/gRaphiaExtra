#' @title testSeuratClusterEnrichment - Test Seurat Clusters for Enrichment
#'
#' @description
#' Performs appropriate tests of enrichment of Seurat clusters for a specified set
#' of genes.
#'
#' @param markers_list Data frame. A data frame containing two columns: a
#'   \code{gene} column and a \code{cluster} column. Typically, output from the
#'   FindAllMarkers() function. If using a custom markers_list object, ensure no
#'   clusters are labelled as "0". If so, rename that cluster to prevent erroneous
#'   enrichment. "0" is a cluster label internally reserved for the background set.
#' @param test_set Data frame or character. Either a data frame including two columns,
#'  "gene" and "group" (e.g., causal genes for different mouse diseases), or a vector of
#'  genes giving the test set (e.g., mouse breast cancer genes).
#' @param background_set Character. A vector of genes that form the background
#'  population for the tests (e.g., all mouse protein coding genes).
#' @param test Character. One of "fisher" (default) for fisher's exact test, "chi.sq"
#'  for the chi squared test, "hyper" for the hypergeometric test, and "conditional"
#'  for partial odds ratios from a multivariable logistic regression. Conditional
#'  is only applicable when groups are included in the test set and there is substantial
#'  overlap between the groups.
#' @param alternative Character. One of "greater" (default), "two.sided", or "less".
#' @param p_adjust Character. Correct for multiple testing using one of "holm",
#'  "hochberg", "hommel", "bonferroni" (default), "BH", "BY", or "fdr" methods.
#'
#' @return A data frame with enrichment test results.
#'
#' @examples
#' # Define Background Set
#' background_set <- c(paste0("gene", 1:20), "gene_causal1", "gene_causal2", "gene_causal3")
#'
#' # Define Markers List
#' # Note: Cluster labels must avoid "0" as specified in your function documentation
#' markers_list <- data.frame(gene = c("gene1", "gene2", "gene_causal1", "gene_causal2",
#'                                     "gene3", "gene4", "gene5", "gene_causal1",
#'                                     "gene6", "gene7", "gene8"),
#'                            cluster = c(rep("Cell Type 1", 4),
#'                                        rep("Cell Type 2", 4),
#'                                        rep("Cell Type 3", 3)),
#'                                        stringsAsFactors = FALSE)
#'
#' # Define Test Set
#' # Data frame with 'gene' and 'group' columns
#' test_set_df <- data.frame(gene = c("gene_causal1", "gene_causal2", "gene_causal3"),
#'                           group = c("Disease_A", "Disease_A", "Disease_B"),
#'                           stringsAsFactors = FALSE)
#'
#' # Example Function Call
#' # One-sided (greater) Fisher's exact test with Bonferroni correction
#' testSeuratClusterEnrichment(markers_list = markers_list,
#'                             test_set = test_set_df,
#'                             background_set = background_set)
#'
#' @export
#' @importFrom dplyr %>%

testSeuratClusterEnrichment <- function(markers_list,
                                        test_set,
                                        background_set,
                                        test = "fisher",
                                        alternative = "greater",
                                        p_adjust = "bonferroni") {
  # global logic
  is_groups <- FALSE

  # check inputs
  checkmate::assert_data_frame(markers_list, min.rows = 1, min.cols = 2, col.names = "unique")
  checkmate::assert_names(colnames(markers_list), must.include = c("gene", "cluster"))
  checkmate::assert_character(markers_list$gene, any.missing = FALSE)
  # make sure clusters are character
  # this is mostly because some users may have already annotated the cluster numbers
  markers_list$cluster <- as.character(markers_list$cluster)
  checkmate::assert_character(markers_list$cluster, any.missing = FALSE)
  checkmate::assert_choice(test, choices = c("fisher", "chi.sq", "hyper", "conditional"))
  if (test == "conditional") {
    checkmate::assert_data_frame(test_set, min.rows = 1, min.cols = 2, col.names = "unique")
    checkmate::assert_names(colnames(test_set), must.include = c("gene", "group"))
    checkmate::assert_character(test_set$gene, any.missing = FALSE)
    checkmate::assert_character(test_set$group, any.missing = FALSE)
  } else {
    checkmate::assert_multi_class(test_set, classes = c("character", "data.frame"))
    if (is.data.frame(test_set)) {
      checkmate::assert_names(colnames(test_set), must.include = c("gene", "group"))
      checkmate::assert_character(test_set$gene, any.missing = FALSE)
      checkmate::assert_character(test_set$group, any.missing = FALSE)
      # update logic
      is_groups <- TRUE
      # clean up dataset for duplicate genes
      test_set <- test_set %>%
        dplyr::distinct(group, gene, .keep_all = TRUE)
    } else {
      checkmate::assert_character(test_set, any.missing = FALSE, min.len = 1, unique = TRUE)
    }
  }
  # get cluster sizes
  # background set must be larger than the largest cluster
  bg_size_check <- table(markers_list$cluster)
  checkmate::assert_character(background_set, any.missing = FALSE, min.len = max(bg_size_check), unique = TRUE)
  if (test == "hyper") {
    checkmate::assert_choice(alternative, choices = c("greater", "less"))
  } else if (test == "fisher") {
    checkmate::assert_choice(alternative, choices = c("greater", "two.sided", "less"))
  } else {
    checkmate::assert_choice(alternative, choices = c("two.sided"), null.ok = TRUE)
  }
  checkmate::assert_choice(p_adjust, choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"))

  # above were checks for input parameters
  # additional checks of input data
  # check inclusion of test set and membership genes in background set
  # number of genes in background set
  background <- length(background_set)
  # membership genes not in background set
  if (any(!markers_list$gene %in% background_set)) {
    message("Removing genes in markers_list not present in the background set.")
    markers_list <- markers_list[markers_list$gene %in% background_set, ]
  }
  # test set genes not in background set
  if (test == "conditional") {
    if (any(!test_set$gene %in% background_set)) {
      message("Removing genes in test set not present in the background set.")
      test_set <- test_set[test_set$gene %in% background_set, ]
    }
  } else {
    if (isTRUE(is_groups)) {
      if (any(!test_set$gene %in% background_set)) {
        message("Removing genes in test set not present in the background set.")
        test_set <- test_set[test_set$gene %in% background_set, ]
      }

      # update user
      message(paste0("Testing ", length(unique(test_set$group)), " groups."))

    } else {
      if (any(!test_set %in% background_set)) {
        message("Removing genes in test set not present in the background set.")
        test_set <- intersect(test_set, background_set)
      }

      # update user
      message(paste0("Testing ", length(test_set), " genes."))
    }
  }

  # set global variables
  group <- NULL
  gene <- NULL
  beta_total <- NULL
  beta_direct <- NULL
  cluster <- NULL
  or <- NULL
  p <- NULL
  par_or <- NULL
  par_p <- NULL
  p_adj <- NULL
  par_p_adj <- NULL

  ### MAIN FUNCTION
  # Initialize a list to store results
  results_list <- list()

  # unique number of clusters
  clusters <- unique(markers_list$cluster)

  # Initialize progress bar
  pb <- utils::txtProgressBar(min = 0, max = length(clusters), style = 3)

  if (test == "conditional") {
    # 2. Build a wide data frame
    # Each column is a binary indicator (1/0)
    model_df <- data.frame(gene = background_set)

    # Add the Outcome (Is it in the specific cluster we are testing?)
    model_df <- dplyr::left_join(model_df, markers_list, by = "gene")
    model_df$cluster[is.na(model_df$cluster)] <- 0

    # Add the Predictors (Is it in each nosology group?)
    for(grp in unique(test_set$group)) {
      model_df[[grp]] <- as.factor(ifelse(model_df$gene %in% test_set$gene[test_set$group == grp], 1, 0))
    }

    # initialise empty lists
    fit <- list()
    fit_simple <- list()

    for (cl_i in seq_along(clusters)) {
      cl <- clusters[cl_i]

      # restructure working df
      temp_df <- model_df
      temp_df$cluster <- ifelse(temp_df$cluster == cl, 1, 0)

      # The formula 'in_cluster ~ .' tells R to use all columns as predictors
      fit[[as.character(cl)]] <- stats::glm(cluster ~ .,
                                            data = temp_df %>% dplyr::select(-gene),               # exclude the 'gene' name column
                                            family = "binomial")

      # Initialize the sub-list for this cluster
      fit_simple[[as.character(cl)]] <- list()

      for (grp in unique(test_set$group)) {
        # 'reformulate' creates: cluster ~ group_name
        safe_grp <- paste0("`", grp, "`")
        form <- stats::as.formula(paste("cluster ~", safe_grp))

        fit_simple[[as.character(cl)]][[as.character(grp)]] <- stats::glm(form,
                                                                          data = temp_df,
                                                                          family = "binomial")
      }

      # Update progress bar
      utils::setTxtProgressBar(pb, cl_i)
    }

    for (cl_i in seq_along(clusters)) {
      cl <- clusters[cl_i]
      # 1. Extract Full Model (Direct Effects)
      stats_full <- as.data.frame(summary(fit[[as.character(cl)]])$coefficients)
      colnames(stats_full) <- c("beta_direct", "StdErr", "z_value", "par_p")
      stats_full$group <- rownames(stats_full)

      # 2. Extract Simple Models (Total Effects)
      # We need to reach back into the logic of how fit_simple was built for this cluster
      simple_list <- list()

      for (grp in unique(test_set$group)) {
        # Re-running or retrieving the simple model for this specific cluster context
        # Note: If you didn't save fit_simple as a nested list, it's best to tidy it here
        s_stats <- summary(fit_simple[[as.character(cl)]][[as.character(grp)]])$coefficients

        simple_list[[as.character(grp)]] <- data.frame(group = grp,
                                                       beta_total = s_stats[2, 1],            # The coefficient for the group
                                                       p = s_stats[2, 4])                     # p-value
      }

      simple_coeffs <- dplyr::bind_rows(simple_list)

      # 3. Merge and Calculate Proportions
      results_list[[as.character(cl)]] <- stats_full %>%
        dplyr::filter(group != "(Intercept)") %>%
        dplyr::mutate(group = gsub("^`|`1$|`$", "", group)) %>%
        dplyr::left_join(simple_coeffs, by = "group") %>%
        dplyr::mutate(cluster = cl,
                      or = exp(beta_total),
                      par_or = exp(beta_direct)) %>%
        dplyr::select(cluster, group, or, p, par_or, par_p)
    }
  } else {
    for (cl_i in seq_along(clusters)) {
      # set update for progress bar
      cl <- clusters[cl_i]

      # conditional on whether groups are in test set
      if (isTRUE(is_groups)) {
        # initialise empty list for results for groups
        results_list_group <- list()

        # genes in cluster
        genes_in_cluster <- unique(markers_list$gene[markers_list$cluster == cl])

        for (grp in unique(test_set$group)) {
          # create contingency table
          # Logic for a, b, c, d
          # genes in group
          test_set_group <- unique(test_set$gene[test_set$group == grp])

          # a: Genes in BOTH cluster AND nosology (the overlap)
          a_genes <- intersect(genes_in_cluster, test_set_group)
          a <- length(a_genes)

          # b: Genes in cluster but NOT in nosology
          b <- length(genes_in_cluster) - a

          # c: Genes in nosology but NOT in this cluster
          c <- length(test_set_group) - a

          # d: Genes in neither (the rest of the universe)
          d <- background - (a + b + c)

          # chi squared test
          if (test == "chi.sq") {
            # Run test
            # Yates correction is applied by default (correct = TRUE)
            test_result <- stats::chisq.test(matrix(c(a, b, c, d), nrow = 2))

            # set to results list
            results_list_group[[as.character(grp)]] <- data.frame(cluster = cl,
                                                                  group = grp,
                                                                  y.cluster_y.test = a,
                                                                  y.cluster_n.test = b,
                                                                  n.cluster_y.test = c,
                                                                  n.cluster_n.test = d,
                                                                  yy_genes = paste(a_genes, collapse = ", "),
                                                                  chi_stat = test_result$statistic,
                                                                  p = test_result$p.value,
                                                                  stringsAsFactors = FALSE)
          }

          # fold enrichment test
          if (test == "hyper") {
            if (alternative == "greater") {
              fe_pval <- stats::phyper(q = a - 1,
                                       m = a + c,
                                       n = background - a - c,
                                       k = a + b,
                                       lower.tail = FALSE)
            } else {
              fe_pval <- stats::phyper(q = a,
                                       m = a + c,
                                       n = background - a - c,
                                       k = a + b,
                                       lower.tail = TRUE)
            }

            # set to results list
            results_list_group[[as.character(grp)]] <- data.frame(cluster = cl,
                                                                  group = grp,
                                                                  sS = a,
                                                                  nS = a + b,
                                                                  sP = a + c,
                                                                  nP = background,
                                                                  sS_genes = paste(a_genes, collapse = ", "),
                                                                  fe = (a/(a + b))/((a + c)/background),
                                                                  p = fe_pval)
          }

          # create matrix
          # test fisher's exact test
          if (test == "fisher") {
            test_result <- stats::fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = alternative)

            # set to results list
            results_list_group[[as.character(grp)]] <- data.frame(cluster = cl,
                                                                  group = grp,
                                                                  y.cluster_y.test = a,
                                                                  y.cluster_n.test = b,
                                                                  n.cluster_y.test = c,
                                                                  n.cluster_n.test = d,
                                                                  yy_genes = paste(a_genes, collapse = ", "),
                                                                  or = test_result$estimate,
                                                                  ci_up = test_result$conf.int[2],
                                                                  ci_lo = test_result$conf.int[1],
                                                                  p = test_result$p.value)
          }
        }

        # bind cluster results
        results_list[[as.character(cl)]] <- dplyr::bind_rows(results_list_group)
      } else {
        # when groups not in test set
        # create contingency table
        # Logic for a, b, c, d
        genes_in_cluster <- unique(markers_list$gene[markers_list$cluster == cl])

        # a: Genes in BOTH cluster AND nosology (the overlap)
        a_genes <- intersect(genes_in_cluster, test_set)
        a <- length(a_genes)

        # b: Genes in cluster but NOT in nosology
        b <- length(genes_in_cluster) - a

        # c: Genes in nosology but NOT in this cluster
        c <- length(test_set) - a

        # d: Genes in neither (the rest of the universe)
        d <- background - (a + b + c)

        # chi squared test
        if (test == "chi.sq") {
          # Run test
          # Yates correction is applied by default (correct = TRUE)
          test_result <- stats::chisq.test(matrix(c(a, b, c, d), nrow = 2))

          # set to results list
          results_list[[as.character(cl)]] <- data.frame(cluster = cl,
                                                         y.cluster_y.test = a,
                                                         y.cluster_n.test = b,
                                                         n.cluster_y.test = c,
                                                         n.cluster_n.test = d,
                                                         yy_genes = paste(a_genes, collapse = ", "),
                                                         chi_stat = test_result$statistic,
                                                         p = test_result$p.value,
                                                         stringsAsFactors = FALSE)
        }

        # fold enrichment test
        if (test == "hyper") {
          if (alternative == "greater") {
            fe_pval <- stats::phyper(q = a - 1,
                                     m = a + c,
                                     n = background - a - c,
                                     k = a + b,
                                     lower.tail = FALSE)
          } else {
            fe_pval <- stats::phyper(q = a,
                                     m = a + c,
                                     n = background - a - c,
                                     k = a + b,
                                     lower.tail = TRUE)
          }

          # set to results list
          results_list[[as.character(cl)]] <- data.frame(cluster = cl,
                                                         sS = a,
                                                         nS = a + b,
                                                         sP = a + c,
                                                         nP = background,
                                                         sS_genes = paste(a_genes, collapse = ", "),
                                                         fe = (a/(a + b))/((a + c)/background),
                                                         p = fe_pval)
        }

        # create matrix
        # test fisher's exact test
        if (test == "fisher") {
          test_result <- stats::fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = alternative)

          # set to results list
          results_list[[as.character(cl)]] <- data.frame(cluster = cl,
                                                         y.cluster_y.test = a,
                                                         y.cluster_n.test = b,
                                                         n.cluster_y.test = c,
                                                         n.cluster_n.test = d,
                                                         yy_genes = paste(a_genes, collapse = ", "),
                                                         or = test_result$estimate,
                                                         ci_up = test_result$conf.int[2],
                                                         ci_lo = test_result$conf.int[1],
                                                         p = test_result$p.value)
        }
      }

      # Update progress bar
      utils::setTxtProgressBar(pb, cl_i)
    }
  }

  # final results
  # bind rows to form results df
  results_df <- dplyr::bind_rows(results_list)
  rownames(results_df) <- NULL

  # correct for multiple testing
  # get adjusted p-values
  results_df$p_adj <- stats::p.adjust(results_df$p, method = p_adjust)
  if (test == "conditional") {
    results_df$par_p_adj <- stats::p.adjust(results_df$par_p, method = p_adjust)

    # rearrange columns to maintain order
    results_df <- results_df %>%
      dplyr::select(cluster, group, or, p, p_adj, par_or, par_p, par_p_adj)
  }

  # close progress bat
  close(pb)

  # return results df
  return(results_df)
}
