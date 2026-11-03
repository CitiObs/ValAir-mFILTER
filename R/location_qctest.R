#' Spatio-Temporal Attributes
#'
#' Checks whether geolocated measurements lie inside inside a region of interest (ROI), whether the implied movement speed is plausible.
#' Assigns a quality flag to each measurement.
#'
#' @details
#'
#' Let \eqn{X_t} denote the measurement at timestamp \eqn{t}, with longitude \eqn{\lambda} and latitude \eqn{\phi}. For each observations, the separation distance \eqn{d_t} (m), the elapsed time \eqn{\Delta t} (s), and instantaneous speed \eqn{v_t} (\eqn{\text{km h}^{-1}})
#' relative to the previous observation \eqn{t-1} are computed as:
#'
#' \deqn{
#' d_t =
#' \begin{cases}
#' \text{D}((\lambda_{t_{-1}}, \phi_{t_{-1}}), (\lambda_t, \phi_t)), & t > t_{-1} \\
#' \text{NA}, & \text{otherwise}
#' \end{cases},
#' \quad
#' \Delta t =
#' \begin{cases}
#' t - t_{-1}, & t > t_{-1} \\
#' \text{NA}, & \text{otherwise}
#' \end{cases},
#' \quad
#' v_t =
#' \begin{cases}
#' \dfrac{d_t}{\Delta t} \times 3.6, & t > t_{-1} \\
#' \text{NA}, & \text{otherwise}
#' \end{cases}
#' }
#'
#' where \eqn{\text{D}(\cdot)} is the geodesic (great-circle) distance on the WGS84 ellipsoid.
#'
#' The quality flag \eqn{f_t^{\text{location}}} is defined as:
#'
#' \deqn{
#' f_t^{\text{location}} =
#' \begin{cases}
#' 0 & \text{Invalid: } (\lambda_t, \phi_t) \notin \text{ROI}, \text{ or } v_t > v_{\text{th}}, \text{ or } d_t = 0, \\
#' 1 & \text{Valid: otherwise}
#' \end{cases}
#' }
#'
#' with:
#' \itemize{
#'   \item \eqn{\text{ROI}}: region defined by `bounding.box`. If `bounding.box`
#'     is `NULL`, the ROI check is skipped.
#'   \item \eqn{v_{\text{th}}}: speed limit (`v.th`).
#' }
#'
#' @param data A data frame containing the dataset to be tested.
#' @param date.column character. Name of the date/time column in `data`.
#' @param id.column character. Name of the column identifying the measurement location ID.
#' @param x.column character. Name of the longitude column in `data`.
#' @param y.column character. Name of the longitude column in `data`.
#' @param v.th numeric. Speed limit in km/h. Values above it are flagged.
#' @param bounding.box numeric. Bounding box of the ROI, given as `c(x.min, x.max, y.min, y.max)` in the same coordinate system as `x.column` and `y.column`. Default is `NULL`.
#'
#' @return The original dataset with:
#' \itemize{
#'   \item \code{flag_location}: integer. Quality flags based on the statistical test. \code{0} = invalid, \code{1} = valid. It is `NA` where the tested value is missing.
#' }
#'
#' @export
location_qctest <- function(data,
                            date.column,
                            id.column = NULL,
                            x.column,
                            y.column,
                            v.th,
                            bounding.box = NULL) {
  # argument checks
  stopifnot(is.character(date.column), date.column %in% colnames(data))
  stopifnot(is.character(x.column), x.column %in% colnames(data))
  stopifnot(is.character(y.column), y.column %in% colnames(data))
  v.th <- ifelse(is.null(v.th), 140, v.th)
  # bounding box filtering
  if (!is.null(bounding.box)) {
    if (length(bounding.box) != 4) {
      stop('The bounding box require xmin. xmax, ymin, and ymxax coordinates!')
    } else {
      data <- data %>%
        filter(.data[[x.column]] >= bounding.box[1], .data[[x.column]] <= bounding.box[2],
               .data[[y.column]] >= bounding.box[3], .data[[y.column]] >= bounding.box[4])
    }
  }
  # location test application
  location_fun <- function(df) {
    df <- df %>%
      arrange(.data[[date.column]])
    df.sf <- st_as_sf(df, coords = c(x.column, y.column), crs = 4326, remove = FALSE)
    df.sf <- df.sf %>%
      mutate(dd = as.numeric(st_distance(geometry, lag(geometry), by_element = TRUE)),
             dt = as.numeric(difftime(.data[[date.column]], lag(.data[[date.column]]), units = "secs")),
             v  = ifelse(dt > 0, 3.6 * dd / dt, NA_real_))
    df <- st_drop_geometry(df.sf)
    #
    df <- df %>%
      mutate(flag_location = case_when(v > v.th | dd == 0 ~ 0,
                                       is.na(dd) & is.na(dt) ~ 1,
                                       v <= v.th ~ 1,
                                       TRUE ~ NA_real_)) %>%
      select(-all_of(c('dd', 'dt', 'v')))
    return(df)
  }
  # Run on groups or entire data
  if (!is.null(id.column)) {
    stopifnot(is.character(id.column), id.column %in% colnames(data))
    data <- data %>%
      group_by(across(all_of(id.column))) %>%
      group_map(~ location_fun(.x), .keep = TRUE) %>%
      bind_rows()
  } else {
    data <- location_fun(data)
  }
  return(as_tibble(data))
}
