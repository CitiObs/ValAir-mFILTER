#' Find Nearby Measurement Points (Nodes and References)
#'
#' Finds the neighboring locations of each node within a given distance, optionally limited to a maximum number of neighbors. The target locations can be moving, i.e. the same identifier can have different coordinates at different timestamps. The neighbors can be static (e.g. reference stations) or moving (e.g. sensor nodes), in which case they are matched to the targets in time.
#'
#' @details
#' For a site \eqn{i} at timestamp \eqn{t}, the selected neighboring sites are:
#' \deqn{
#' \mathcal{J}^{\ast}_t =
#' \left\{
#' j \in \mathcal{J}
#' \;\middle|\;
#' d_{ijt} \le \rho
#' \;\wedge\;
#' \mathrm{rank}(d_{ijt}) \leq n_{\text{maxneigh}}
#' \right\}
#' }
#'
#' with:
#' \itemize{
#'   \item \eqn{\mathcal{J}}: set of all candidate neighboring locations in `rdata`,
#'   \item \eqn{\mathcal{J}^{\ast}_t}: subset of \eqn{\mathcal{J}} that satisfies the spatial and, for moving neighbors, temporal criteria,
#'   \item \eqn{d_{ijt}}: geographic distance in m between target \eqn{i} at timestamp \eqn{t} and neighbor \eqn{j} at timestamp \eqn{t^{\ast}_j},
#'   \item \eqn{\rho} (`maxradii`): radius of influence, i.e. the maximum distance at which a site is considered a neighbor,
#'   \item \eqn{n_{\text{maxneigh}}} (`maxneigh`): maximum number of neighbors to keep for each target.
#' }
#'
#' Neighbors are ranked in ascending order of \eqn{d_{ijt}}, so the closest ones are kept. If `maxneigh = NULL`, all locations within \eqn{\rho} are retained. If `exclude.colocated = TRUE`, locations at exactly the same position as the target (\eqn{d_{ijt} = 0}) are excluded.\cr
#'
#' The timestamp \eqn{t^{\ast}_j} of the neighbor depends on its type (`static`) and on the matching mode (`match`):
#' \itemize{
#'   \item `static = TRUE`: the neighbors have fixed coordinates and are available at all timestamps, so \eqn{t^{\ast}_j = t}.
#'   \item `static = FALSE`: the neighbors can change location with time and are matched to the target in time according to `match`:
#'   \itemize{
#'     \item `match = "exact"`: only neighbors observed at the same timestamp as the target are considered, so \eqn{t^{\ast}_j = t}.
#'     \item `match = "near"`: the observation of neighbor \eqn{j} with the closest timestamp to \eqn{t} is used, provided that it lies within the time tolerance:
#'     \deqn{
#'     |t - t^{\ast}_j| \le \tau,
#'     }
#'     where \eqn{\tau} (`time.tolerance`) is the maximum allowed time difference, in minutes
#'   }
#' }
#'
#' If `aggregation = TRUE` (only available for static neighbors, `static = TRUE`), the values of the selected neighbors are averaged for each target at its timestamp \eqn{t} using the modified Shepard weights:
#'
#' \deqn{
#' w_{ij} = \left( \frac{\rho - d_{ijt}}{\rho \, d_{ijt}} \right)^{2}
#' \qquad
#' \hat{z}_i(t) = \frac{\sum_{j \in \mathcal{J}^{\ast}_t} w_{ij} \, z_j(t)}{\sum_{j \in \mathcal{J}^{\ast}_t} w_{ij}}
#' }
#'
#' with:
#' \itemize{
#'   \item \eqn{z_j(t)}: value of the variable at neighbor \eqn{j} at timestamp \eqn{t},
#'   \item \eqn{\mathcal{J}^{\ast}_t}: in this formula, the subset of selected neighbors with a non-missing value at time \eqn{t}.
#' }
#'
#' Neighbors with \eqn{d_{ijt} = 0} have an infinite weight, so if any is present, the average of their values is returned.
#'
#' @param data A data frame containing the dataset to be tested.
#' @param date.column character. Name of the date/time column in `data` and `rdata`. In `rdata`, it is required if `static = FALSE` or `aggregation = TRUE`.
#' @param id.column character. Name of the column identifying the measurement location in `data`.
#' @param rdata A data frame containing additional locations (e.g. reference stations or validated nodes) used as neighbors.
#' @param rdata.column character. Name of the column(s) in `rdata` to aggregate. Only used if `aggregation = TRUE`. If `NULL`, all numeric columns of `rdata` are aggregated. Default is `NULL`.
#' @param rid.column character. Name of the column identifying the measurement location in `rdata`. If it is the same as `id.column`, the neighbor identifier is renamed with the prefix `r_` in the output.
#' @param x.column character. Name of the longitude column in `data` and `rdata`.
#' @param y.column character. Name of the latitude column in `data` and `rdata`.
#' @param static logical. If `TRUE`, the locations in `rdata` are fixed (e.g. reference stations). If `FALSE`, they can change with time and are matched to the targets in time according to `match`. Default is `TRUE`.
#' @param match character. Matching mode, used only if `static = FALSE`: `"exact"` if `data` and `rdata` have the same temporal resolution, or `"near"` to match the nearest timestamp within `time.tolerance`. Default is `"exact"`.
#' @param time.tolerance numeric. Maximum allowed time difference in minutes Used only if `static = FALSE` and `match = "near"`. Default is `0`.
#' @param maxradii numeric. Radius of influence in km. Default is `30`.
#' @param maxneigh numeric. Maximum number of neighbors to keep for each target and timestamp. Default is `NULL`.
#' @param exclude.colocated logical. If `TRUE`, neighbors at exactly the same position as the target (\eqn{d_{ijt} = 0}) are excluded. Default is `FALSE`.
#' @param aggregation logical. If `TRUE`, the data of the selected neighbors are spatially aggregated at the timestamps of `data`. It is only available for static neighbors (`static = TRUE`), i.e. reference stations. Default is `FALSE`.
#'
#' @return A list with two data frames, or with both elements `NULL` if there are no neighbors:
#' \itemize{
#'   \item \code{data}: if `aggregation = FALSE`, the subset of `rdata` that belongs to the selected neighbors (all their rows if `static = TRUE`, and only the matched observations if `static = FALSE`). If `aggregation = TRUE`, the aggregated data of the selected neighbors for each target, date and position.
#'   \item \code{neighbors}: spatial information of the selected neighbors for each target: the identifier, date, longitude and latitude of the target, followed by the identifier (`rid.column`), date, longitude and latitude (prefixed with `r_`) of the neighbor and the separation distance \code{dij} in m.
#' }
#'
#' @seealso \code{\link{spatial_similarity_qctest}} and \code{\link{chain_similarity_qctest}} use this helper function.
#'
#' @keywords internal
find_neighbors <- function(data,
                           date.column,
                           id.column,
                           rdata,
                           rdata.column = NULL,
                           rid.column,
                           x.column,
                           y.column,
                           static = TRUE,
                           match = c('exact', 'near'),
                           time.tolerance = 0,
                           maxradii = 30,
                           maxneigh = NULL,
                           exclude.colocated = FALSE,
                           aggregation = FALSE) {
  # argument checks
  match <- match.arg(match)
  stopifnot(is.character(id.column), id.column %in% colnames(data))
  stopifnot(is.character(date.column), date.column %in% colnames(data))
  stopifnot(is.character(rid.column), rid.column %in% colnames(rdata))
  stopifnot(is.character(x.column), x.column %in% colnames(data), x.column %in% colnames(rdata))
  stopifnot(is.character(y.column), y.column %in% colnames(data), y.column %in% colnames(rdata))
  stopifnot(is.logical(static), length(static) == 1, !is.na(static))
  stopifnot(is.numeric(maxradii), length(maxradii) == 1, maxradii > 0)
  stopifnot(is.logical(exclude.colocated), length(exclude.colocated) == 1, !is.na(exclude.colocated))
  stopifnot(is.logical(aggregation), length(aggregation) == 1, !is.na(aggregation))
  if (!is.null(maxneigh)) {
    stopifnot(is.numeric(maxneigh), length(maxneigh) == 1, maxneigh >= 1)
  }
  if (!static) {
    stopifnot(date.column %in% colnames(rdata))
    stopifnot(is.numeric(time.tolerance), length(time.tolerance) == 1, !is.na(time.tolerance), time.tolerance >= 0)
  }
  if (aggregation) {
    stopifnot("`aggregation = TRUE` is only available for static neighbors (`static = TRUE`)" = static)
    stopifnot(date.column %in% colnames(rdata))
    # by default, all numeric columns of rdata are aggregated
    if (is.null(rdata.column)) {
      rdata.column <- rdata %>%
        select(where(is.numeric)) %>%
        colnames() %>%
        setdiff(c(rid.column, x.column, y.column, date.column))
    }
    stopifnot(is.character(rdata.column), rdata.column %in% colnames(rdata), length(rdata.column) > 0)
  }
  #
  if (id.column == rid.column) {
    rdata <- rdata %>%
      rename(!!paste0("r_", rid.column) := all_of(rid.column))
    rid.column <- paste0("r_", rid.column)
  }
  # maximum radii in m and time tolerance in min converted to s (an exact match allows no time difference)
  maxradii <- maxradii * 1000
  tolerance <- if (match == 'exact') 0 else time.tolerance * 60
  # target locations (moving: one row per id, date and position)
  info.target <- data %>%
    select(all_of(c(id.column, date.column, x.column, y.column))) %>%
    filter(!is.na(.data[[x.column]]), !is.na(.data[[y.column]])) %>%
    distinct()
  # neighbors locations (static: one row per id and position; moving: one row per id, date and position)
  info.neigh <- rdata %>%
    select(all_of(c(rid.column, if (!static) date.column, x.column, y.column))) %>%
    filter(!is.na(.data[[x.column]]), !is.na(.data[[y.column]])) %>%
    distinct() %>%
    rename(!!paste0("r_", x.column) := all_of(x.column),
           !!paste0("r_", y.column) := all_of(y.column))
  if (!static) {
    info.neigh <- info.neigh %>%
      rename(!!paste0("r_", date.column) := all_of(date.column))
  }
  if (nrow(info.neigh) == 0 || nrow(info.target) == 0) {
    return(list(data = NULL, neighbors = NULL))
  }
  if (static) {
    # separation distance of the static neighbors (rows = neighbors, columns = targets)
    sp.target <- st_as_sf(info.target, crs = 4326, coords = c(x.column, y.column))
    sp.neigh <- st_as_sf(info.neigh, crs = 4326, coords = c(paste0("r_", x.column), paste0("r_", y.column)))
    distances <- matrix(as.numeric(sf::st_distance(sp.neigh, sp.target)), nrow = nrow(info.neigh))
    # candidate pairs (target, neighbor) within the maximum radii
    hit <- which(distances <= maxradii, arr.ind = TRUE)
    info.neigh <- bind_cols(info.target[hit[, "col"], , drop = FALSE],
                            info.neigh[hit[, "row"], , drop = FALSE]) %>%
      mutate(!!paste0("r_", date.column) := .data[[date.column]],
             dij = distances[hit])
  } else {
    # moving neighbors: observation times in seconds, sorted for each neighbor
    to_sec <- function(x) as.numeric(as.POSIXct(x))
    # index of the closest value of the sorted vector y for each value of x
    nearest_index <- function(x, y) {
      i <- findInterval(x, y)
      i0 <- pmax(i, 1)
      i1 <- pmin(i + 1, length(y))
      ifelse(abs(x - y[i0]) <= abs(y[i1] - x), i0, i1)
    }
    neigh.list <- info.neigh %>%
      arrange(.data[[paste0("r_", date.column)]]) %>%
      group_by(across(all_of(rid.column))) %>%
      group_split()
    neigh.sec <- lapply(neigh.list, function(x) to_sec(x[[paste0("r_", date.column)]]))
    t.target <- to_sec(info.target[[date.column]])
    # closest observation in time of each neighbor for all targets, within the time tolerance
    info.neigh <- lapply(seq_along(neigh.list), function(k) {
      idx <- nearest_index(t.target, neigh.sec[[k]])
      ok <- abs(neigh.sec[[k]][idx] - t.target) <= tolerance
      if (!any(ok)) {
        return(NULL)
      }
      bind_cols(info.target[ok, , drop = FALSE], neigh.list[[k]][idx[ok], , drop = FALSE])
    }) %>%
      bind_rows()
    if (nrow(info.neigh) == 0) {
      return(list(data = NULL, neighbors = NULL))
    }
    # separation distance between each target and its matched neighbor observation
    sp.target <- st_as_sf(info.neigh, crs = 4326, coords = c(x.column, y.column))
    sp.neigh <- st_as_sf(info.neigh, crs = 4326, coords = c(paste0("r_", x.column), paste0("r_", y.column)))
    info.neigh <- info.neigh %>%
      mutate(dij = as.vector(sf::st_distance(sp.neigh, sp.target, by_element = TRUE))) %>%
      filter(.data[["dij"]] <= maxradii)
  }
  # selection of neighbors for each target (id, date, position), closest first
  if (exclude.colocated) {
    info.neigh <- info.neigh %>%
      filter(.data[["dij"]] > 0)
  }
  info.neigh <- info.neigh %>%
    group_by(across(all_of(c(id.column, date.column, x.column, y.column)))) %>%
    arrange(.data[["dij"]], .by_group = TRUE)
  if (!is.null(maxneigh)) {
    info.neigh <- info.neigh %>%
      slice_head(n = maxneigh)
  }
  info.neigh <- info.neigh %>%
    ungroup() %>%
    select(all_of(c(id.column, date.column, x.column, y.column, rid.column,
                    paste0("r_", date.column), paste0("r_", x.column), paste0("r_", y.column), "dij")))
  if (nrow(info.neigh) == 0) {
    return(list(data = NULL, neighbors = NULL))
  }
  # aggregation of neighbors data with Shepard weights (static neighbors only)
  if (aggregation) {
    # weighted mean ignoring missing values
    wmean <- function(v, d) {
      ok <- !is.na(v)
      if (!any(ok)) return(NA_real_)
      v <- v[ok]
      d <- d[ok]
      # co-located neighbors have an infinite weight
      if (any(d == 0)) return(mean(v[d == 0]))
      w <- ((maxradii - d) / (maxradii * d))^2
      sum(v * w) / sum(w)
    }
    # neighbors values at the timestamps of the data
    df.values <- rdata %>%
      select(all_of(c(rid.column, date.column, rdata.column))) %>%
      distinct()
    df.neigh <- info.neigh %>%
      inner_join(df.values, by = c(rid.column, date.column)) %>%
      group_by(across(all_of(c(id.column, date.column, x.column, y.column)))) %>%
      summarise(across(all_of(rdata.column), ~ wmean(.x, .data[["dij"]])),
                .groups = "drop")
  } else if (static) {
    # subset rdata to selected neighbors
    df.neigh <- rdata %>%
      filter(.data[[rid.column]] %in% info.neigh[[rid.column]])
  } else {
    # subset rdata to the matched neighbor observations
    df.neigh <- rdata %>%
      semi_join(info.neigh %>%
                  select(all_of(c(rid.column, paste0("r_", date.column)))) %>%
                  distinct() %>%
                  rename(!!date.column := all_of(paste0("r_", date.column))),
                by = c(rid.column, date.column))
  }
  #
  return(list(data = df.neigh, neighbors = info.neigh))
}
