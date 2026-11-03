#' Spatial Correlation Test
#'
#' Evaluates whether a measurement location is correlated with its neighboring locations using Pearson's R correlation coefficient. Assigns a quality flag and optionally creates a column of valid values.
#'
#' @details
#' Measurements with sub-hourly resolution are aggregated to hourly intervals using the median operator, applied to both the measured values and the location coordinates. If the data are already hourly, the aggregation step is skipped.\cr
#'
#' Let \eqn{X_{i\tau}} denote the original sub-hourly measurements at node \eqn{i} and timestamp \eqn{\tau}, and \eqn{(\lambda_{i\tau}, \phi_{i\tau})} the location of node \eqn{i} at that timestamp. The hourly series \eqn{X_{it}} and the hourly location \eqn{(\lambda_{it}, \phi_{it})} are:
#'
#' \deqn{
#' X_{it} = \operatorname{median} \left\{ X_{i\tau} \mid \tau \in h \right\}, \quad
#' \lambda_{it} = \operatorname{median} \left\{ \lambda_{i\tau} \mid \tau \in h \right\}, \quad
#' \phi_{it} = \operatorname{median} \left\{ \phi_{i\tau} \mid \tau \in h \right\}
#' }
#'
#' where \eqn{t} is an hourly interval and \eqn{\tau \in h} are the sub-hourly
#' timestamps within that hour.
#'
#' The test can be applied to:
#' \itemize{
#'   \item the entire dataset (`historical = FALSE`), or
#'   \item consecutive rolling windows \eqn{\mathcal{T}} of `ndays` days (`historical = TRUE`). One window is evaluated per day. Each window ends at 23:00 of day \eqn{n} and contains the current day and the preceding `ndays - 1` days:
#'     \deqn{
#'     \mathcal{T} = \left[ t_n - (T - 1)\,\Delta t,\; t_n \right],
#'     }
#'     where \eqn{t_n} is 23:00 of day \eqn{n}, \eqn{\Delta t} is the temporal resolution (1 hour) and \eqn{T = n_{\text{days}} \cdot 24} is the number of hourly time steps in the window. The flag obtained for a window is assigned to all timestamps of day \eqn{n}.
#' }
#'
#' For each node \eqn{i} and temporal window, the separation distance \eqn{d_{ijt}} between \eqn{i} and every other site \eqn{j} is computed at each hour \eqn{t}.
#' Site \eqn{j} is a neighbor of \eqn{i} at hour \eqn{t} only if \eqn{d_{ijt} \le \rho}, where \eqn{\rho} (`corr.radii`) is the radius of influence.
#' Since the locations change with time, the neighbors of \eqn{i} can differ from hour to hour.\cr
#'
#' The pairwise Pearson correlation coefficient \eqn{R_{ij}} between \eqn{i} and \eqn{j} is then computed using only the hours at which both sites are neighbors:
#'
#' \deqn{
#' R_{ij}(\mathcal{T}) =
#' \frac{
#'   \sum\limits_{t \in \mathcal{T}^{\ast}}
#'   (X_{it} - \bar{X}_i^{\ast})(X_{jt} - \bar{X}_j^{\ast})
#' }{
#'   \sqrt{\sum\limits_{t \in \mathcal{T}^{\ast}} (X_{it} - \bar{X}_i^{\ast})^2}
#'   \sqrt{\sum\limits_{t \in \mathcal{T}^{\ast}} (X_{jt} - \bar{X}_j^{\ast})^2}
#' },
#' }
#'
#' where:
#' \itemize{
#'   \item \eqn{\mathcal{T}^{\ast} \subseteq \mathcal{T}} denotes the timestamps for which both \eqn{X_{it}} and \eqn{X_{jt}} are available (non-missing) and the sites are within the radius of influence (\eqn{d_{ijt} \le \rho}):
#'   \deqn{
#'   \mathcal{T}^{\ast} = \left\{ t \in \mathcal{T} \mid X_{it}, X_{jt} \text{ available} \;\land\; d_{ijt} \le \rho \right\}
#'   }
#'   \item \eqn{\bar{X}_i^{\ast}} and \eqn{\bar{X}_j^{\ast}} are the temporal
#'   means of \eqn{X_{it}} and \eqn{X_{jt}} over \eqn{\mathcal{T}^{\ast}}:
#'   \deqn{
#'   \bar{X}_i^{\ast} = \frac{1}{N_{\mathcal{T}^{\ast}}}
#'   \sum_{t \in \mathcal{T}^{\ast}} X_{it}, \quad
#'   \bar{X}_j^{\ast} = \frac{1}{N_{\mathcal{T}^{\ast}}}
#'   \sum_{t \in \mathcal{T}^{\ast}} X_{jt}
#'   }
#'   \item \eqn{N_{\mathcal{T}^{\ast}}}: number of valid paired observations
#'   within \eqn{\mathcal{T}^{\ast}}.
#' }
#'
#' \strong{Correlation threshold.} The correlation \eqn{R_{ij}} is compared with a threshold \eqn{R_{\text{th},j}} that depends on `corr.method`:
#' \itemize{
#'   \item `"standard"`: the same fixed threshold is used for all neighbors, \eqn{R_{\text{th},j} = R_{\text{th}}} (`corr.th`).
#'   \item `"lut"`: the threshold is the \emph{expected correlation} of neighbor \eqn{j}, obtained from a look-up function \eqn{R_{exp}(s_t, d_{ijt})} that depends on the season \eqn{s_t} of hour \eqn{t}, and the separation distance \eqn{d_{ijt}} (in km). See \code{\link{expected_correlation}}.
#'   The expected correlation is evaluated at every hour at which \eqn{j} is within the radius of influence, and the threshold of neighbor \eqn{j} is the median over those hours:
#'   \deqn{
#'   R_{\text{th},j} = \operatorname{median} \left\{ R_{exp}(s_t, d_{ijt}) \mid t \in \mathcal{T}, \; d_{ijt} \le \rho \right\}
#'   }
#' }
#'
#' Two counts are computed:
#'
#' \deqn{
#' N_j(\rho) = \left| \left\{ j \mid N_{\mathcal{T}^{\ast}} \ge n_{\text{min}} \right\} \right|,
#' }
#' \deqn{
#' N_{R_{ij} > R_{\text{th}}}(\rho) = \left| \left\{ j \mid R_{ij}(\mathcal{T}) > R_{\text{th},j} \;\land\; N_{\mathcal{T}^{\ast}} \ge n_{\text{min}} \right\} \right|,
#' }
#'
#' where:
#' \itemize{
#'   \item \eqn{N_j(\rho)}: number of sites that are neighbors of site \eqn{i} (within the radius \eqn{\rho}) for at least \eqn{n_{\text{min}}} paired hours, i.e. enough data to compute a reliable correlation. If it is `0`, the site is isolated.
#'   \item \eqn{N_{R_{ij} > R_{\text{th}}}(\rho)}: number of those neighbors whose correlation with site \eqn{i} exceeds their threshold \eqn{R_{\text{th},j}}. If it is `0`, the site is not correlated with any neighbor. A neighbor whose threshold cannot be computed (`NA`) is not counted as correlated.
#'   \item \eqn{R_{\text{th},j}}: correlation threshold for neighbor \eqn{j}; constant (`corr.th`) for `corr.method = "standard"` and neighbor-specific (expected correlation) for `corr.method = "lut"`.
#'   \item \eqn{n_{\text{min}}}: minimum number of paired hourly observations required to compute a reliable correlation and to run the test. It is calculated as \eqn{n_{\text{min}} = p \times n_{\text{days}} \times 24}, where \eqn{p} (`compl.corr.duration`) is the required data completeness, expressed as a fraction of the hours in the window, and \eqn{n_{\text{days}}} (`ndays`) is the window length in days.
#'   \item \eqn{d_{ijt}}: distance between sites \eqn{i} and \eqn{j} at hour \eqn{t}.
#'   \item \eqn{\rho} (`corr.radii`): radius of influence. Only sites within this distance of site \eqn{i} at hour \eqn{t} (\eqn{d_{ijt} \le \rho}) are considered neighbors at that hour.
#' }
#'
#' The quality flag \eqn{f_i^{\text{correlation}}} for site \eqn{i} is defined as:
#'
#' \deqn{
#' f_i^{\text{correlation}} =
#' \begin{cases}
#' 0 & \text{Insufficient data: } N_i < n_{\text{min}}, \\
#' 1 & \text{Isolated site: } N_j(\rho) < 1, \\
#' 2 & \text{No correlation with neighbors: } N_{R_{ij} > R_{\text{th}}} < 1, \\
#' 3 & \text{Correlation with neighbors: } N_{R_{ij} > R_{\text{th}}} \ge 1.
#' \end{cases}
#' }
#'
#' where \eqn{N_i} is the number of valid hourly values of site \eqn{i} within the window.
#'
#' When the test is applied dynamically (e.g., using consecutive rolling windows), the same quality flag applies to all timestamps within day \eqn{n}:
#'
#' \deqn{
#' f_{it}^{\text{correlation}}(n) = f_i^{\text{correlation}}
#' }
#'
#' Otherwise, a single quality flag is assigned to all timestamps \eqn{t} and native timestamps \eqn{\tau \in t}:
#'
#' \deqn{
#' \begin{array}{cc}
#' f_{it}^{\text{correlation}} = f_i^{\text{correlation}}, &
#' f_{i\tau}^{\text{correlation}} = f_{it}^{\text{correlation}}
#' \end{array}
#' }
#'
#' @param data A data frame containing the data to be tested.
#' @param date.column character. Name of the date/time column in `data` and `rdata`.
#' @param data.column character. Name of the column in `data` to test.
#' @param id.column character. Name of the column identifying the measurement location in `data`.
#' @param rdata A data frame containing additional locations (e.g. reference stations) used as neighbors.
#' @param rdata.column character. Name of the column in `rdata` used in the test.
#' @param rid.column character. Name of the column identifying the measurement location in `rdata`.
#' @param x.column character. Name of the longitude column in `data` and `rdata`.
#' @param y.column character. Name of the latitude column in `data` and `rdata`.
#' @param compl.corr.duration numeric. Minimum data completeness required within the temporal window to run the test, expressed as a fraction of the hours in the window. The number of hours in the window is `ndays * 24`.Default is `0.75`.
#' @param corr.radii numeric. Radius of influence in km. Two locations are neighbors at a given hour if their separation distance at that hour is not larger than this value. Default is `10`.
#' @param corr.method character. Method used to define the correlation threshold. One of:
#' \itemize{
#'   \item `"standard"`: a fixed threshold, `corr.th`, is used for all neighbors.
#'   \item `"lut"`: a neighbor-specific threshold is used, given by the expected correlation (median over the hours within the radius) from `expected_correlation()` as a function of season and separation distance. `corr.th` is ignored.
#' }
#' Default is `"standard"`.
#' @param corr.th numeric. Correlation threshold above which a neighbor counts as correlated. Used only if `corr.method = "standard"`. Default is `0.75`.
#' @param historical logical. If `TRUE`, the test is applied in consecutive rolling windows; otherwise, it is applied to the entire `data`. Default is `TRUE`.
#' @param ndays numeric. Window length in days if `historical = TRUE`. It is also used to compute the minimum number of observations required (`compl.corr.duration * ndays * 24`). Default is `30`.
#' @param new logical. If `TRUE`, adds a column `<param>_correlation` containing only the valid values; values with a flag other than `3` are replaced with `NA`. Default is `TRUE`.
#' @param param character. Base name of the new column (e.g., `PM2.5` gives `PM2.5_correlation`). Default is `"PM2.5"`. Used only if `new = TRUE`.
#'
#' @return The original dataset with:
#' \itemize{
#'   \item \code{flag_correlation}: numeric. Quality flags based on the statistical test. \code{0} = insufficient data, \code{1} = isolated site, \code{2} = no correlation with neighbors, \code{3} = correlated with neighbors. It is `NA` where the tested value is missing.
#'   \item \code{<param>_correlation}: numeric. Valid values only (if `new = TRUE`).
#' }
#'
#' @export
spatial_correlation_qctest <- function(data,
                                       date.column,
                                       data.column,
                                       id.column,
                                       rdata = NULL,
                                       rdata.column = NULL,
                                       rid.column = NULL,
                                       x.column,
                                       y.column,
                                       compl.corr.duration= 0.75,
                                       corr.radii = 10,
                                       corr.method = 'standard',
                                       corr.th = 0.75,
                                       historical = TRUE,
                                       ndays = 30,
                                       new = TRUE,
                                       param = 'PM2.5') {
  # Argument checks
  stopifnot(all(sapply(c(date.column, data.column, id.column, x.column, y.column), is.character)))
  stopifnot(all(c(date.column, data.column, id.column, x.column, y.column) %in% colnames(data)))
  stopifnot(is.numeric(compl.corr.duration), compl.corr.duration >= 0, compl.corr.duration <= 1)
  stopifnot(is.numeric(corr.radii), !is.null(corr.radii))
  stopifnot(is.logical(historical))
  #
  data <- data %>%
    mutate(!!id.column := as.character(.data[[id.column]]))
  #
  # correlation radii to m
  corr.radii <- corr.radii * 1000
  n.min <- compl.corr.duration*ndays*24
  # correlation test
  data %>%
    group_by(across(all_of(id.column))) %>%
    group_map(~ {
      id.target <- .y[[1]]
      df <- .x
      #
      correlation_process_fun <- function(df.target) {
        if (sum(!is.na(df.target[[data.column]])) < n.min) return(0)
        timestamps.range <- range(df.target[[date.column]])
        # neighbor data and incorporate rdata if exists
        if (!is.null(rdata)) {
          # arguments for rdata
          stopifnot(all(sapply(list(rdata.column, date.column, rid.column), function(x) !is.null(x) && is.character(x))))
          stopifnot(all(c(rdata.column, date.column, rid.column) %in% colnames(rdata)))
          #
          rdata <- rdata %>%
            select(all_of(c(date.column, rid.column, rdata.column, x.column, y.column))) %>%
            rename(!!id.column := !!rid.column,
                   !!data.column := !!rdata.column) %>%
            mutate(!!id.column := match.fun(paste0("as.", class(data[[id.column]])[1]))(.data[[id.column]])) %>%
            select(all_of(c(date.column, id.column, x.column, y.column, data.column)))
          df.neigh <- bind_rows(data %>% filter(.data[[id.column]] != id.target) , rdata) %>% as_tibble()
        } else {
          df.neigh <- data %>%
            filter(.data[[id.column]] != id.target)
        }
        df.neigh <- df.neigh %>%
          filter(.data[[date.column]] >= timestamps.range[1],
                 .data[[date.column]] <= timestamps.range[2])
        #
        # hourly aggregation of the data and of the location
        df.target.agg <- df.target %>%
          mutate(!!date.column := floor_date(.data[[date.column]], 'hour')) %>%
          group_by(across(all_of(c(id.column, date.column)))) %>%
          summarise(!!x.column := median(.data[[x.column]], na.rm = TRUE),
                    !!y.column := median(.data[[y.column]], na.rm = TRUE),
                    !!data.column := median(.data[[data.column]], na.rm = TRUE), .groups = "drop") %>%
          arrange(.data[[date.column]]) %>%
          group_by(.data[[id.column]]) %>%
          complete(!!date.column := seq(min(.data[[date.column]], na.rm=TRUE), max(.data[[date.column]], na.rm=TRUE), by = 'hour')) %>%
          ungroup()
        df.neigh.agg <- df.neigh %>%
          mutate(!!date.column := floor_date(.data[[date.column]], 'hour')) %>%
          group_by(across(all_of(c(id.column, date.column)))) %>%
          summarise(!!x.column := median(.data[[x.column]], na.rm = TRUE),
                    !!y.column := median(.data[[y.column]], na.rm = TRUE),
                    !!data.column := median(.data[[data.column]], na.rm = TRUE), .groups = "drop") %>%
          arrange(.data[[date.column]]) %>%
          group_by(.data[[id.column]]) %>%
          complete(!!date.column := seq(min(.data[[date.column]], na.rm=TRUE), max(.data[[date.column]], na.rm=TRUE), by = 'hour')) %>%
          ungroup()
        #
        # omit neighbor data outside the temporal window
        df.neigh.agg <- df.neigh.agg %>%
          group_by(.data[[id.column]]) %>%
          group_split() %>%
          lapply(function(x) {
            if (sum(!is.na(x[[data.column]])) < n.min) NULL else x
          }) %>%
          rlist::list.clean(fun = is.null) %>%
          bind_rows()
        #
        if (nrow(df.neigh.agg) == 0) return(1)
        #
        # target location in each hour
        info.target <- df.target.agg %>%
          filter(!is.na(.data[[x.column]]), !is.na(.data[[y.column]])) %>%
          select(all_of(c(date.column, x.column, y.column))) %>%
          rename(!!paste0(x.column, '_target') := all_of(x.column),
                 !!paste0(y.column, '_target') := all_of(y.column))
        #
        # neighbor location in the same hour as the target
        df.neigh <- df.neigh.agg %>%
          filter(!is.na(.data[[x.column]]), !is.na(.data[[y.column]]), !is.na(.data[[data.column]])) %>%
          inner_join(info.target, by = date.column)
        if (nrow(df.neigh) == 0) return(1)
        sp.target <- st_as_sf(df.neigh %>% select(all_of(c(paste0(x.column, '_target'), paste0(y.column, '_target')))), coords = c(paste0(x.column, '_target'), paste0(y.column, '_target')), crs = 4326L)
        sp.neigh <- sf::st_as_sf(df.neigh %>% select(all_of(c(x.column, y.column))), coords = c(x.column, y.column), crs = 4326L)
        #
        # calculation of separation distance in each hour and selection of neighbors
        distances <- as.vector(round(st_distance(sp.target, sp.neigh, by_element = TRUE), 1))
        df.neigh <- df.neigh %>%
          mutate(dij = distances) %>%
          filter(dij <= corr.radii)
        #
        if (nrow(df.neigh) == 0) return(1)
        #
        # neighboring and target data in matrix format (neighbor values only in the hours where it is within the radius of the target)
        mat.neigh <- pivot_wider(df.neigh %>% select(all_of(c(date.column, id.column, data.column))),
                                 names_from = all_of(id.column),
                                 values_from = all_of(data.column))
        mat.neigh <- left_join(df.target.agg %>% select(all_of(c(date.column, data.column))), mat.neigh, by = date.column) %>%
          select(-all_of(c(date.column, data.column)))
        #
        # Compute pairwise correlation between target and neighbors
        mat.target <- replicate(ncol(mat.neigh), df.target.agg[[data.column]])
        mat.completeness <- apply(mat.target + mat.neigh, 2, function(x) sum(!is.na(x)))
        completeness.stat <- tibble(!!id.column := names(mat.completeness), completeness := as.vector(mat.completeness))
        pairwise.correlation <- mapply(cor, as.data.frame(mat.neigh), as.data.frame(mat.target), method='pearson', use='pairwise.complete.obs')
        correlation.stat <- tibble(!!id.column := names(pairwise.correlation), correlation := as.vector(pairwise.correlation))
        df.correlation <- left_join(correlation.stat, completeness.stat, by=id.column)
        #
        # Comparisons
        if (corr.method == 'standard') {
          # fixed threshold for all neighbors
          df.comparison <- df.correlation %>%
            filter(completeness >= n.min)
          if (nrow(df.comparison) == 0) return(1)
          df.comparison <- df.comparison %>%
            filter(correlation > corr.th)
        } else {
          # neighbor-specific expected correlation (season + median distance)
          df.expected.correlation <- df.neigh %>%
            mutate(Season = quarter(.data[[date.column]], fiscal_start = 12)) %>%
            mutate(expected_correlation = expected_correlation(season = Season, distance = dij/1000, relax.factor = 0.25)) %>%
            select(-Season)
          df.expected.correlation <- df.expected.correlation %>%
            group_by(node) %>%
            summarise(expected_correlation = median(expected_correlation, na.rm = TRUE))
          #
          df.comparison <- left_join(df.correlation,
                                     df.expected.correlation,
                                     by = id.column) %>%
            filter(completeness >= n.min)
          if (nrow(df.comparison) == 0) return(1)
          df.comparison <- df.comparison %>%
            filter(correlation > expected_correlation)
        }
        df.comparison <- df.correlation %>%
          filter(completeness >= n.min)
        if (nrow(df.comparison) == 0) return(1)
        df.comparison <- df.comparison %>%
          filter(correlation > corr.th)
        if (nrow(df.comparison) == 0) return(2)
        return(3)
      }
      #
      # check for (not) historical data
      if (historical) {
        stopifnot(is.numeric(ndays), ndays > 0, !is.null(ndays))
        # run in temporal windows
        days.intervals <- ymd_hms(paste(unique(as_date(df[[date.column]])), "23:00:00"))
        days <- as_date(floor_date(days.intervals, unit = 'day'))
        correlation.flags <- runner::runner(x = df,
                                            f = function(x) correlation_process_fun(x),
                                            k = paste(ndays, "days"),
                                            idx = df[[date.column]],
                                            at = days.intervals,
                                            simplify = TRUE)
        daily.correlation.flag <- tibble(Date = days, flag_correlation = correlation.flags)
        df <- df %>%
          mutate(Date = as_date(.data[[date.column]])) %>%
          left_join(daily.correlation.flag, by = "Date") %>%
          select(-Date)
        df <- df %>%
          mutate(flag_correlation = ifelse(is.na(.data[[data.column]]), NA_real_, flag_correlation))
      } else {
        correlation.flags <- correlation_process_fun(df)
        df <- df %>%
          mutate(flag_correlation = ifelse(is.na(.data[[data.column]]), NA_real_, correlation.flags))
      }
      if (new == TRUE) {
        df <- df %>%
          mutate(!!paste0(param, "_correlation") := ifelse(flag_correlation != 3 | is.na(flag_correlation), NA_real_, .data[[data.column]]))
      }
      return(as_tibble(df))
    }, .keep = TRUE) %>%
    bind_rows()
}
