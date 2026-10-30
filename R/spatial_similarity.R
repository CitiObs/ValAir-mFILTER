#' Spatial Similarity Function
#'
#' Evaluates whether measurements are consistent with those of neighbouring higher-accuracy nodes and assigns a quality flag.
#'
#' @details
#' For each time step \eqn{t_i}, the temporal window \eqn{\mathcal{T}} is defined as:
#' \deqn{
#' \mathcal{T} = [t_i - D, t_i],
#' }
#' where \eqn{D} (`similarity.duration`) is the duration of the window, in hours.
#'
#' Let \eqn{X_{it}} be the measurements at site \eqn{i} and \eqn{Y_{jt}} those of a neighbouring higher-accuracy node \eqn{j}. The similarity indicator \eqn{S_{ij}(\mathcal{T})} is computed with one of
#' two methods (`method`):
#'
#' \itemize{
#'   \item Median Absolute Deviation (`method = 'mad'`):
#'   \deqn{
#'     MAD_{ij}(\mathcal{T}^{\ast}) = \operatorname{median}_{t \in \mathcal{T}^{\ast}}
#'     \left\{
#'       \left| r_{ijt} - \operatorname{median}_{t \in \mathcal{T}^{\ast}} r_{ijt} \right|
#'     \right\},
#'     \quad
#'     r_{ijt} = \left| X_{it} - Y_{jt} \right|
#'   }
#'
#'   \item Normalized Root Mean Square Deviation (`method = 'nrmsd'`):
#'   \deqn{
#'     nRMSD_{ij}(\mathcal{T}^{\ast}) =
#'     \frac{
#'       \sqrt{\frac{1}{n_{\mathcal{T}^{\ast}}} \sum_{t \in \mathcal{T}^{\ast}} (X_{it} - Y_{jt})^2}
#'     }{
#'       \frac{1}{n_{\mathcal{T}^{\ast}}} \sum_{t \in \mathcal{T}^{\ast}} Y_{jt}
#'     }
#'   }
#' }
#'
#' where \eqn{\mathcal{T}^{\ast} \subseteq \mathcal{T}} denotes the timestamps for which both \eqn{X_{it}} and \eqn{Y_{jt}} are available (non-missing), and \eqn{n_{\mathcal{T}^{\ast}}} is their number.
#'
#' \eqn{S_{ij}(\mathcal{T})} is compared against a threshold
#' \eqn{S_{\text{th}}} (`similarity.threshold`). Two counts are then derived
#' for each window:
#' \deqn{
#' N_{j} = \left| \left\{ j : n_{\mathcal{T}^{\ast}} \ge n_{\text{min}} \right\} \right|,
#' }
#' \deqn{
#' N_{S_{ij} < S_{\text{th}}} = \left| \left\{ j : S_{ij}(\mathcal{T}) < S_{\text{th}} \;\land\; n_{\mathcal{T}^{\ast}} \ge n_{\text{min}} \right\} \right|,
#' }
#'
#' where:
#' \itemize{
#'   \item \eqn{N_{j}}: number of higher-accuracy nodes with enough paired data.
#'   \item \eqn{N_{S_{ij} < S_{\text{th}}}}: number of those nodes whose similarity is below the threshold.
#'   \item \eqn{S_{\text{th}}} (`similarity.threshold`): similarity threshold.
#'   \item \eqn{n_{\text{min}}}: minimum number of paired hourly observations required in the window. It is calculated as \eqn{n_{\text{min}} = p \times D}, where \eqn{p} (`compl.similarity.duration`) is the required data completeness, expressed as a fraction of the hours in the window, and \eqn{D} (`similarity.duration`) is the window length in hours.
#' }
#'
#' The quality flag \eqn{f_{it}^{\text{similarity}}} for site \eqn{i} at each timestamp \eqn{t} is defined as:
#'
#' \deqn{
#' f_{it}^{\text{similarity}} = \begin{cases}
#' 0 & \text{Insufficient data: } N_{j} < 1, \\
#' 1 & \text{Isolated site: no higher-accuracy nodes available}, \\
#' 2 & \text{No similarity with nearby higher-accuracy nodes: } N_{S_{ij} < S_{\text{th}}} < 1, \\
#' 3 & \text{Similarity with nearby higher-accuracy nodes: } N_{S_{ij} < S_{\text{th}}} \geq 1.
#' \end{cases}
#' }
#'
#' @param data A data frame containing the data to be tested.
#' @param date.column character. Name of the date/time column in `data` and `rdata`.
#' @param data.column character. Name of the column in `data` to test.
#' @param id.column character. Name of the column identifying the measurement location in `data`.
#' @param rdata A data frame containing the data of the higher-accuracy (reference) stations.
#' @param rdata.column character. Name of the column in `rdata` used in the test.
#' @param rid.column character. Name of the column identifying the measurement location in `rdata`.
#' @param similarity.duration numeric. Length of the temporal window, in hours. Default is `72`.
#' @param compl.similarity.duration numeric. Minimum data completeness required within `similarity.duration` to run the test, expressed as a fraction of the hours in the window (greater than 0 and at most 1). Default is `0.75`.
#' @param method character. Similarity indicator, one of `"mad"` and `"nrmsd"`. Default is `"mad"`.
#' @param similarity.threshold numeric. Threshold of similarity indicator. Defaults are `5` for `method = "mad"` and `0.68` for `method = "nrmsd"`.
#'
#' @return The input data with an additional column, \code{flag_similarity}, containing the quality flag: \code{0} = insufficient data,
#' \code{1} = isolated site, \code{2} = no similarity with nearby nodes, \code{3} = similarity with nearby nodes. It is `NA` where the tested value is missing.
#'
#' @export
spatial_similarity <- function(data,
                               date.column,
                               data.column,
                               id.column = NULL,
                               rdata,
                               rdata.column,
                               rid.column,
                               similarity.duration = 72,
                               compl.similarity.duration = 0.75,
                               method = 'mad',
                               similarity.threshold = NULL) {
  # Early exit if rdata is NULL or not available
  if (is.null(rdata)) {
    data <- data %>% mutate(flag_similarity = 1)
    return(data)
  }
  # Argument checks
  stopifnot(is.character(date.column), date.column %in% colnames(data))
  stopifnot(is.character(data.column), data.column %in% colnames(data))
  stopifnot(is.character(rdata.column), rdata.column %in% colnames(rdata))
  stopifnot(is.character(rid.column), rid.column %in% colnames(rdata))
  stopifnot(is.numeric(similarity.duration), similarity.duration > 0)
  stopifnot(is.numeric(compl.similarity.duration),
            compl.similarity.duration > 0, compl.similarity.duration <= 1)
  stopifnot(method %in% c('mad', 'nrmsd'))
  if (!is.null(id.column)) {
    stopifnot(is.character(id.column), id.column %in% colnames(data))
  }
  if (is.null(similarity.threshold)) similarity.threshold <- if (method == "mad") 5 else 0.68
  stopifnot(is.numeric(similarity.threshold))
  #
  # minimum number of paired hourly observations in the window: n_min = p * D
  n.min <- ceiling(compl.similarity.duration * similarity.duration)
  #
  ind.nona <- which(!is.na(data[[data.column]]))
  if (length(ind.nona) == 0) {
    data <- data %>% mutate(flag_similarity = NA_real_)
  } else {
    start.nona <- data[[date.column]][ind.nona[1]]
    end.nona <- data[[date.column]][ind.nona[length(ind.nona)]]
    data.sub <- data %>%
      filter(.data[[date.column]] >= start.nona, .data[[date.column]] <= end.nona)
    if (sum(!is.na(data.sub[[data.column]])) < n.min) {
      data <- data %>% mutate(flag_similarity = 0)
    } else {
      data.neigh <- rdata %>%
        filter(.data[[date.column]] >= start.nona, .data[[date.column]] <= end.nona)
      #
      mat.neigh <- tidyr::pivot_wider(data.neigh %>% select(all_of(c(date.column, rid.column, rdata.column))) %>% distinct(),
                                      names_from = all_of(rid.column),
                                      values_from = all_of(rdata.column))
      mat.target.neigh <- left_join(data.sub %>% select(all_of(c(date.column, data.column))), mat.neigh, by = date.column)
      mat.neigh <- mat.target.neigh %>% select(-all_of(c(date.column, data.column)))
      mat.target <- replicate(ncol(mat.neigh), mat.target.neigh[[data.column]])
      #
      if (method == 'mad') {
        mat.similarity <- as_tibble(abs(mat.target - mat.neigh))
        similarity.stat <- apply(mat.similarity, 2, function(x) {
          slider::slide_index_dbl(x, data.sub[[date.column]],
                                  ~ median(abs(.x - median(.x, na.rm = TRUE)), na.rm = TRUE),
                                  .before = lubridate::hours(similarity.duration), .after = 0, .complete = FALSE)
        })
      } else if (method == 'nrmsd') {
        mat.similarity <- as_tibble(mat.target - mat.neigh)
        roll.rmsd <- apply(mat.similarity, 2, function(x) {
          slider::slide_index_dbl(x, data.sub[[date.column]],
                                  ~ sqrt(mean(.x^2, na.rm = TRUE)),
                                  .before = lubridate::hours(similarity.duration), .after = 0, .complete = FALSE)
        })
        roll.mean <- apply(mat.neigh, 2, function(x) {
          slider::slide_index_dbl(x, data.sub[[date.column]],
                                  ~ mean(.x, na.rm = TRUE),
                                  .before = lubridate::hours(similarity.duration), .after = 0, .complete = FALSE)
        })
        similarity.stat <- roll.rmsd / roll.mean
      }
      # number of valid paired observations in each window
      validity.mat <- as_tibble(!is.na(mat.target + mat.neigh))
      completeness.stat <- apply(validity.mat, 2, function(x) {
        slider::slide_index_dbl(x, data.sub[[date.column]],
                                ~ sum(.x),
                                .before = lubridate::hours(similarity.duration), .after = 0, .complete = FALSE)
      })
      #
      mat.comparison <- bind_cols(!!date.column := data.sub[[date.column]], as_tibble(similarity.stat))
      mat.comparison <- tidyr::pivot_longer(mat.comparison,
                                            cols = -all_of(c(date.column)),
                                            names_to = rid.column,
                                            values_to = 'Statistic')
      compl.comparison <- bind_cols(!!date.column := data.sub[[date.column]], as_tibble(completeness.stat))
      compl.comparison <- tidyr::pivot_longer(compl.comparison,
                                              cols = -all_of(c(date.column)),
                                              names_to = rid.column,
                                              values_to = 'Completeness')
      mat.comparison <- left_join(mat.comparison, compl.comparison, by = c(date.column, rid.column))
      mat.comparison.results <- mat.comparison %>%
        group_by(across(all_of(date.column))) %>%
        summarise(n_id_available = n(),
                  any_available = any(!is.na(Completeness) & Completeness >= n.min),
                  any_similar = any(!is.na(Completeness) & !is.na(Statistic) & Completeness >= n.min & Statistic < similarity.threshold),
                  flag_similarity = case_when(any_similar ~ 3,
                                              !any_similar & any_available ~ 2,
                                              n_id_available < 1 ~ 1,
                                              TRUE ~ 0), .groups = "drop")
      data.sub <- left_join(data.sub, mat.comparison.results %>% select(all_of(c(date.column, 'flag_similarity'))), by = date.column) %>%
        select(all_of(c(date.column, 'flag_similarity')))
      data <- left_join(data, data.sub, by = date.column)
    }
  }
  return(data)
}
