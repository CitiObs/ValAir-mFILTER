#' Spatial Similarity Test
#'
#' Compares measurements at a target location with those of nearby reference stations, assigns a quality flag, and optionally creates a column of valid values.
#'
#' @details
#' This function is a \strong{high-level wrapper} that:
#' \itemize{
#'   \item identifies nearby reference stations using \code{\link{find_neighbors}} function,
#'   \item computes the spatial similarity statistic and the quality flags using \code{\link{spatial_similarity}} function, and
#'   \item maps the resulting flags back to the native temporal resolution.
#' }
#'
#' Measurements with sub-hourly resolution are aggregated to hourly values using the median operator. For data that are already hourly, the aggregation leaves the values unchanged.
#'
#' Let \eqn{X_{i\tau}} denote the original sub-hourly measurements at node \eqn{i} and timestamp \eqn{\tau}. The hourly series \eqn{X_{it}} is:
#'
#' \deqn{
#' X_{it} = \operatorname{median}\{ X_{i\tau} \mid \tau \in t \},
#' }
#' where \eqn{t} is an hourly interval.
#'
#' The test is applied on the hourly time series within rolling windows \eqn{\mathcal{T} = [t_i - D, t_i]}, where \eqn{D} (`similarity.duration`) is the window length in hours. For each window, the similarity indicator (see \code{\link{spatial_similarity}}) is computed between site
#' \eqn{i} and each neighboring reference station \eqn{j}. A site \eqn{i} is considered similar to a reference station if it:
#'
#' \itemize{
#'   \item is within distance \eqn{\rho} (`similarity.radii`),
#'   \item has sufficient paired data in the window, i.e. at least \eqn{n_{\text{min}} = p \times D} paired hourly values, where \eqn{p} (`compl.similarity.duration`) is the required data completeness, expressed as a fraction of the hours in the window,
#'   \item has a similarity statistic below the threshold (`similarity.threshold`).
#' }
#'
#' The quality flag \eqn{f_{it}^{\text{similarity}}} derived for each hour \eqn{t} is uniformly applied to all native timestamps \eqn{\tau \in t}:
#'
#' \deqn{
#' f_{i\tau}^{\text{similarity}} = f_{it}^{\text{similarity}}.
#' }
#'
#' @param data A data frame containing the data to be tested.
#' @param date.column character. Name of the date/time column in `data` and `rdata`.
#' @param data.column character. Name of the column in `data` to test.
#' @param id.column character. Name of the column identifying the measurement location in `data`.
#' @param rdata A data frame containing the data of the reference stations.
#' @param rdata.column character. Name of the column in `rdata` used in the test.
#' @param rid.column character. Name of the column identifying the measurement location in `rdata`.
#' @param x.column character. Name of the longitude column in `data` and `rdata`.
#' @param y.column character. Name of the latitude column in `data` and `rdata`.
#' @param similarity.radii numeric. Radius of influence in km. Default is `30`.
#' @param similarity.duration numeric. Length of the rolling window, in hours. Default is `72`.
#' @param compl.similarity.duration numeric. Minimum data completeness required within `similarity.duration` to run the test, expressed as a fraction of the hours in the window. Default is `0.75`.
#' @param method character. One of `"mad"` (default) or `"nrmsd"`. See \code{\link{spatial_similarity}}. Default is `"mad"`.
#' @param similarity.threshold numeric. The upper threshold of the calculated statistical index. Defaults are `5` for `method = "mad"` and `0.68` for `method = "nrmsd"`. See \code{\link{spatial_similarity}}.
#' @param new logical. If `TRUE`, adds a column `<param>_similarity` containing only the valid values; values with a flag other than `3` are replaced with `NA`. Default is `FALSE`.
#' @param param character. Base name of the new column (e.g., `PM2.5` gives `PM2.5_similarity`). Default is `"PM2.5"`. Used only if `new = TRUE`.
#'
#' @return The original dataset with:
#' \itemize{
#'   \item \code{flag_similarity}: numeric. Quality flags based on the statistical test. \code{0} = insufficient data, \code{1} = isolated site, \code{2} = non-similar with nearby reference stations, \code{3} = similar with nearby reference stations. It is `NA` where the tested value is missing.
#'   \item \code{<param>_similarity}: numeric. Valid values only (if `new = TRUE`).
#' }
#'
#' @export
spatial_similarity_qctest <- function(data,
                                      date.column,
                                      data.column,
                                      id.column = NULL,
                                      rdata,
                                      rdata.column,
                                      rid.column,
                                      x.column,
                                      y.column,
                                      similarity.radii = 30,
                                      similarity.duration = 168,
                                      compl.similarity.duration = 0.75,
                                      method = 'mad',
                                      similarity.threshold = NULL,
                                      new = FALSE,
                                      param = 'PM2.5') {
  # arguments check
  method <- match.arg(method)
  stopifnot(is.character(data.column), length(data.column) == 1, data.column %in% colnames(data))
  stopifnot(is.logical(new), length(new) == 1, !is.na(new))
  if (new) {
    if (is.null(param)) param <- data.column
    stopifnot(is.character(param), length(param) == 1, !is.na(param))
  }
  # Early exit if rdata is NULL
  if (is.null(rdata)) {
    data <- data %>%
      mutate(flag_similarity = NA_real_)
    if (new) data <- data %>% mutate(!!paste0(param, "_similarity") := NA_real_)
    return(data)
  }
  # arguments check
  stopifnot(is.character(date.column), date.column %in% colnames(data))
  stopifnot(is.character(rdata.column), rdata.column %in% colnames(rdata))
  stopifnot(is.character(rid.column), rid.column %in% colnames(rdata))
  stopifnot(is.character(x.column), x.column %in% colnames(data), x.column %in% colnames(rdata))
  stopifnot(is.character(y.column), y.column %in% colnames(data), y.column %in% colnames(rdata))
  stopifnot(is.numeric(similarity.radii), similarity.radii > 0)
  stopifnot(is.numeric(similarity.duration), similarity.duration > 0)
  stopifnot(is.numeric(compl.similarity.duration), compl.similarity.duration > 0, compl.similarity.duration <= 1)
  if (!is.null(id.column)) {
    stopifnot(is.character(id.column), id.column %in% colnames(data))
  }
  if (is.null(similarity.threshold)) {
    similarity.threshold <- if (method == 'mad') 5 else 0.68
  }
  stopifnot(is.numeric(similarity.threshold), length(similarity.threshold) == 1)
  #
  data <- data %>%
    mutate(!!id.column := as.character(.data[[id.column]]))
  rdata <- rdata %>%
    mutate(!!rid.column := as.character(.data[[rid.column]]))
  #
  similarity_fun <- function(df) {
    # temporal aggregation of the data and rdata
    df.agg <- df %>%
      mutate(!!date.column := lubridate::floor_date(.data[[date.column]], 'hour')) %>%
      group_by(across(all_of(c(id.column, date.column)))) %>%
      summarise(!!data.column := median(.data[[data.column]], na.rm = TRUE),
                !!x.column := median(.data[[x.column]], na.rm = TRUE),
                !!y.column := median(.data[[y.column]], na.rm = TRUE),
                .groups = "drop") %>%
      arrange(.data[[date.column]])
    rdata.agg <- rdata %>%
      mutate(!!date.column := lubridate::floor_date(.data[[date.column]], "hour")) %>%
      group_by(across(all_of(c(rid.column, x.column, y.column, date.column)))) %>%
      summarise(!!rdata.column := median(.data[[rdata.column]], na.rm = TRUE), .groups = "drop")
    #
    # find nearest neighbors (static references), aggregated at the timestamps of the target
    neigh.info <- find_neighbors(data = df.agg,
                                 date.column = date.column,
                                 id.column = id.column,
                                 rdata = rdata.agg,
                                 rdata.column = rdata.column,
                                 rid.column = rid.column,
                                 x.column = x.column,
                                 y.column = y.column,
                                 static = TRUE,
                                 maxradii = similarity.radii,
                                 aggregation = TRUE)
    # no neighbors available
    if (is.null(neigh.info$data)) {
      df <- df %>% mutate(flag_similarity = 1)
    } else {
      # spatial similarity with the aggregated neighbors (a single series)
      rdata.agg <- neigh.info$data %>%
        mutate(!!rid.column := 1)
      #
      spatial.similarity.results <- spatial_similarity(data = df.agg,
                                                       date.column = date.column,
                                                       data.column = data.column,
                                                       id.column = id.column,
                                                       rdata = rdata.agg,
                                                       rdata.column = rdata.column,
                                                       rid.column = rid.column,
                                                       similarity.duration = similarity.duration,
                                                       compl.similarity.duration = compl.similarity.duration,
                                                       method = method,
                                                       similarity.threshold = similarity.threshold) %>%
        select(all_of(c(date.column, "flag_similarity")))
      df <- df %>%
        mutate(timestamp_h = lubridate::floor_date(.data[[date.column]], "hour")) %>%
        left_join(spatial.similarity.results, by = c("timestamp_h" = date.column)) %>%
        select(-timestamp_h) %>%
        mutate(flag_similarity = ifelse(is.na(.data[[data.column]]), NA_real_, flag_similarity))
    }
    if (new) {
      df <- df %>%
        mutate(!!paste0(param, "_similarity") := ifelse(.data[['flag_similarity']] != 3 | is.na(.data[['flag_similarity']]), NA_real_, .data[[data.column]]))
    }
    return(as_tibble(df))
  }
  if (!is.null(id.column)) {
    data <- data %>% group_by(across(all_of(id.column))) %>%
      group_map(~ similarity_fun(.x), .keep = TRUE) %>%
      bind_rows()
  } else {
    data <- similarity_fun(data)
  }
  return(as_tibble(data))
}
