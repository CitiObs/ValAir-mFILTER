#' Constant Value Test
#'
#' Detects flatline behavior, i.e. consecutive measurements that remain almost constant within a temporal window, and assigns a quality flag. Optionally creates a column of valid values.
#'
#' @details
#' Let \eqn{X_t} denote the measurements at a location. For each time step
#' \eqn{t_i}, the temporal window \eqn{\mathcal{T}} is defined as:
#' \deqn{
#' \mathcal{T} = [t_i - D, t_i],
#' }
#' where \eqn{D} (`persist.duration`) is the duration of the window, in minutes.
#'
#' Missing values are excluded from the calculation. Let \eqn{\mathcal{T}^{\ast} \subseteq \mathcal{T}} be the set of available
#' (non-missing) values in the window, and \eqn{n_{\mathcal{T}^{\ast}}} their number. The variability within the window is measured by the range:
#' \deqn{
#' \text{Range}_{\mathcal{T}^{\ast}} =
#' \max_{t \in \mathcal{T}^{\ast}}(X_t) - \min_{t \in \mathcal{T}^{\ast}}(X_t).
#' }
#'
#' The quality flag \eqn{f_{t}^{\text{constant}}} is then defined as:
#' \deqn{
#' f_{t}^{\text{constant}} =
#' \begin{cases}
#' 0 & \text{Insufficient data: } n_{\mathcal{T}^{\ast}} < n_{min}, \\
#' 1 & \text{Constant: } n_{\mathcal{T}^{\ast}} \geq n_{min} \text{ and } \text{Range}_{\mathcal{T}^{\ast}} \leq \text{Range}_{min}, \\
#' 2 & \text{Non-constant: } n_{\mathcal{T}^{\ast}} \geq n_{min} \text{ and } \text{Range}_{\mathcal{T}^{\ast}} > \text{Range}_{min},
#' \end{cases}
#' }
#' \itemize{
#'   \item \eqn{n_{min}}: minimum number of available (non-missing) values required to evaluate the test. It is calculated from the data completeness fraction \eqn{p} (`compl.persist.duration`), the window duration \eqn{D} (`persist.duration`) and the temporal resolution \eqn{\Delta t} of \eqn{X_t}, as \eqn{n_{min} = p D / \Delta t}.
#'   \item \eqn{\text{Range}_{min}} (`min.variation`): minimum acceptable range within the window.
#' }
#'
#' @param data A data frame containing the dataset to be tested
#' @param date.column character. Name of the date/time column in `data`.
#' @param data.column character. Name of the column in `data` to test.
#' @param id.column character. Name of the column identifying the measurement location ID. If `NULL`, the test is applied to the whole dataset; otherwise, it is applied to each individual location. Default is `NULL`.
#' @param persist.duration numeric. Length of the rolling window, in minutes.
#' @param compl.persist.duration numeric. Minimum data completeness required within `persist.duration` to run the test, expressed as a fraction. Default is `0.75`.
#' @param min.variation numeric. Minimum range below which values are considered constant. Default is `0.1`.
#' @param new logical. If TRUE, adds a column `<param>_constant_value` containing the valid values; flagged values are replaced with with `NA`. Default is `FALSE`.
#' @param param character. Base name of the new column (e.g. `PM2.5` gives `PM2.5_constant_value`). Default is `"PM2.5"`. Used only if `new = TRUE`.
#'
#' @return The original dataset with:
#' \itemize{
#'   \item \code{flag_constant_value}: integer. Quality flags based on the statistical test, \code{0} = insufficient data, \code{1} = constant, \code{2} = non-constant. It is `NA` where the tested value is missing.
#'   \item \code{<param>_constant_value}: numeric. Valid values only (if `new = TRUE`).
#' }
#'
#' @export
constant_value_qctest <- function(data,
                                  date.column,
                                  data.column,
                                  id.column = NULL,
                                  persist.duration,
                                  compl.persist.duration,
                                  min.variation,
                                  new = TRUE,
                                  param = 'PM2.5') {
  # Argument checks
  stopifnot(is.character(date.column), date.column %in% colnames(data))
  stopifnot(is.character(data.column), length(data.column) == 1, data.column %in% colnames(data))
  stopifnot(persist.duration > 0, !is.null(persist.duration))
  stopifnot(compl.persist.duration >= 0, compl.persist.duration <= 1, !is.null(compl.persist.duration))
  stopifnot(min.variation >= 0)
  # detect constant value periods
  constant_value_fun <- function(df) {
    df <- df %>%
      arrange(.data[[date.column]])
    resolution <- get_mode(round(diff(as.numeric(df[[date.column]])))[round(diff(as.numeric(df[[date.column]]))) > 0])
    n.min <- ceiling(compl.persist.duration*floor(persist.duration*60/resolution) + 1)
    df <- df  %>%
      mutate(roll_range = slider::slide_index_dbl(.data[[data.column]], .data[[date.column]],function(x) {
        if (all(is.na(x))) {
          NA_real_
        } else {
          max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
        }
      }, .before = lubridate::minutes(persist.duration), .after = 0,.complete = FALSE),
      n = slider::slide_index_dbl(.data[[data.column]], .data[[date.column]], function(x) {
        sum(!is.na(x))
      }, .before = lubridate::minutes(persist.duration), .after = 0,.complete = FALSE),
      flag_constant_value = case_when(n < n.min ~ 0,
                                      roll_range <= min.variation & n >= n.min ~ 1,
                                      roll_range > min.variation & n >= n.min ~ 2,
                                      TRUE ~ NA_real_)) %>%
      select(-all_of(c('roll_range', 'n')))
    if (new) {
      df <- df %>%
        mutate(!!paste0(param, "_constant_value") := ifelse(.data[["flag_constant_value"]] %in% c(0, 1) | is.na(.data[["flag_constant_value"]]), NA_real_, .data[[data.column]]))
    }
    return(df)
  }
  # Run on groups or entire data
  if (!is.null(id.column)) {
    stopifnot(is.character(id.column), id.column %in% colnames(data))
    data <- data %>%
      group_by(across(all_of(id.column))) %>%
      group_map(~ constant_value_fun(.x), .keep = TRUE) %>%
      bind_rows()
  } else {
    data <- constant_value_fun(data)
  }
  return(as_tibble(data))
}
