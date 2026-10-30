#' Chain Similarity Test
#'
#' Checks whether measurements at a location are similar to those of neighboring higher-accuracy nodes (official reference stations or already
#' validated nodes), assigns a quality flag, and optionally creates a column of valid values.
#'
#' @details
#' The test evaluates whether measurements at node \eqn{i} are consistent with nearby higher-quality nodes through a \strong{two-stage similarity chain}. Its purpose is to
#' extend spatial validation beyond reference stations: nodes that were already validated against a reference can act as secondary comparison anchors.
#'
#' If node \eqn{i} has a reference station \eqn{j} nearby, it is compared directly with it (Step 1). Otherwise, it can still be evaluated indirectly through a node \eqn{i_{\ast}} that has already been
#' validated against a reference (Step 2).
#'
#' The similarity indicator \eqn{S} is computed with \code{\link{spatial_similarity}} using the method defined by `method` (`"nrmsd"` or `"mad"`) and is compared with the
#' threshold \eqn{S_{\text{th}}} (`similarity.threshold`). Each step has its own radius of influence (`similarity.radii`), window duration (`similarity.duration`), data completeness (`compl.similarity.duration`), method (`method`) and threshold (`similarity.threshold`).
#'
#' \strong{Step 1: Similarity with the nearest reference stations}
#'
#' For each node \eqn{i}, the nearby reference stations \eqn{j} are identified within a radius of influence \eqn{\rho_1} (`similarity.radii[1]`). Their time series are
#' linearly interpolated to the timestamps of node \eqn{i}. The similarity is then computed with \code{\link{spatial_similarity}} using a rolling window
#' \eqn{\mathcal{T}_1 = [t_i - D_1, t_i]}{T1 = [t_i - D1, t_i]} of duration \eqn{D_1}{D1} hours (`similarity.duration[1]`) and a minimum number of paired observations
#' \eqn{n_{\min,1} = p_1 \times D_1}{n_min1 = p1 * D1}, where \eqn{p_1}{p1} (`compl.similarity.duration[1]`) is the required data completeness, expressed as a fraction of the hours in the window.
#'
#' Measurements of node \eqn{i} that are similar to the nearby reference stations (\eqn{S_{ij} < S_{\text{th},1}}{S_ij < S_th1}) receive flag \code{4} and the node becomes a \strong{high-confidence node} \eqn{i_{\ast}},
#' serving as a comparison point in Step 2. Only measurements of node \eqn{i} with
#' \eqn{f_{it} \in \{0,1\}} in Step 1 (insufficient data, or no reference station \eqn{j} within \eqn{\rho_1}) proceed to Step 2.
#'
#' \strong{Step 2: Similarity with validated nodes}
#'
#' The measurements of the high-confidence nodes \eqn{i_{\ast}} are temporally matched with the measurements of \eqn{i}.
#' For each timestamp \eqn{t_i} of node \eqn{i} and each validated node \eqn{i_{\ast}}, the observation of \eqn{i_{\ast}} closest in time to \eqn{t_i} is selected if:
#' \deqn{
#' |t_{i_{\ast}} - t_i| \le \Delta t,
#' }{|t_i* - t_i| <= Delta t,}
#' with \eqn{\Delta t} the time tolerance (`time.tolerance`).
#' If \eqn{i_{\ast}} has no observation within \eqn{[t_i - \Delta t, t_i + \Delta t]}{[t_i - Delta t, t_i + Delta t]}, it is not used as a neighbor of \eqn{i} at \eqn{t_i}.
#' Each validated node contributes at most one observation per \eqn{t_i}, while several validated nodes can be matched to the same \eqn{t_i}.
#'
#' Only the validated nodes within the radius \eqn{\rho_2} (`similarity.radii[2]`) at the matched timestamp are kept. Their values are averaged at each \eqn{t_i},
#' so that \eqn{i} is compared with a single neighbor time series:
#' \deqn{
#' \bar{X}_{\ast}(t_i) = \frac{1}{|\mathcal{I}^{\ast}_{t_i}|} \sum_{i_{\ast} \in \mathcal{I}^{\ast}_{t_i}} X_{i_{\ast}}(t_{i_{\ast}}),
#' }{z*(t_i) = mean of z_i*(t_i*) over the validated nodes matched to t_i,}
#' where \eqn{\mathcal{I}^{\ast}_{t_i}}{I*_t_i} is the set of validated nodes matched to \eqn{t_i} with a non-missing value, \eqn{|\mathcal{I}^{\ast}_{t_i}|}{|I*_t_i|} their number,
#' and \eqn{X_{i_{\ast}}(t_{i_{\ast}})}{z_i*(t_i*)} the value of node \eqn{i_{\ast}} at its matched observation.
#'
#' The similarity between \eqn{i} and \eqn{\bar{z}_{\ast}}{z*} is then evaluated as in Step 1, using a rolling window
#' \eqn{\mathcal{T}_2 = [t_i - D_2, t_i]}{T2 = [t_i - D2, t_i]} of duration \eqn{D_2}{D2} hours (`similarity.duration[2]`) and a minimum number of paired observations
#' \eqn{n_{\min,2} = p_2 \times D_2}{n_min2 = p2 * D2}, where \eqn{p_2}{p2} is `compl.similarity.duration[2]`.
#'
#' \strong{Chain similarity flag}
#'
#' The flag \eqn{f_{it}^{\text{chain similarity}}} combines both steps. The comparison with a reference station (Step 1) has priority.
#' The comparison with neighboring validated nodes (Step 2) is used only for measurements without a valid Step 1 result, i.e., no reference station
#' within \eqn{\rho_1} or insufficient data (\eqn{f_{it} \in \{0,1\}}).
#'
#' \deqn{
#' f_{it}^{\text{chain similarity}} =
#' \begin{cases}
#' 0 & \text{Insufficient data: }
#'     N_j < 1 \text{ if } \exists j,
#'     \text{ or } N_{i_{\ast}} < 1
#'     \text{ if } \nexists j \land \exists i_{\ast}, \\
#' 1 & \text{Isolated site: }
#'     \nexists j \land \nexists i_{\ast}, \\
#' 2 & \text{No similarity (direct or indirect evaluation): }
#'     S_{ij}(\mathcal{T}_1^{\ast}) \ge S_{\text{th},1} \text{ if } \exists j,
#'     \text{ or } S_{i i_{\ast}}(\mathcal{T}_2^{\ast}) \ge S_{\text{th},2}
#'     \text{ if } \nexists j \land \exists i_{\ast}, \\
#' 3 & \text{Indirect evaluation (validated node) with similarity: }
#'     S_{i i_{\ast}}(\mathcal{T}_2^{\ast}) < S_{\text{th},2}, \\
#' 4 & \text{Direct evaluation (reference station) with similarity: }
#'     S_{ij}(\mathcal{T}_1^{\ast}) < S_{\text{th},1}.
#' \end{cases}
#' }{f = 0 (insufficient data), 1 (isolated site), 2 (no similarity), 3 (similar to a validated node), 4 (similar to a reference station)}
#'
#' where:
#' \itemize{
#'   \item \eqn{\exists j}: at least one reference station within \eqn{\rho_1}.
#'   \item \eqn{\exists i_{\ast}}: at least one validated node within \eqn{\rho_2} with an observation within \eqn{\Delta t} of \eqn{t_i}.
#'   \item \eqn{N_j}, \eqn{N_{i_{\ast}}}: number of reference stations and of validated nodes with at least \eqn{n_{\min,1}}{n_min1} and \eqn{n_{\min,2}}{n_min2} paired observations, respectively.
#'   \item \eqn{S_{\text{th},1}}{S_th1}, \eqn{S_{\text{th},2}}{S_th2} (`similarity.threshold`): thresholds of the similarity indicator in Step 1 and Step 2. If `NULL`, each is `0.68` for `"nrmsd"` and `5` for `"mad"`.
#' }
#'
#' The flag is assigned at each native timestamp \eqn{t_i} of node \eqn{i}, since no temporal aggregation is applied.
#'
#' @param data A data frame containing the data to be tested.
#' @param date.column character. Name of the date/time column in `data` and `rdata`.
#' @param data.column character. Name of the column in `data` to test.
#' @param id.column character. Name of the column identifying the measurement location in `data`.
#' @param rdata A data frame containing the data of the reference stations. If `NULL`, the flags (and the new column, if `new = TRUE`) are set to `NA`.
#' @param rdata.column character. Name of the column in `rdata` used in the test.
#' @param rid.column character. Name of the column identifying the measurement location in `rdata`.
#' @param x.column character. Name of the longitude column in `data` and `rdata`.
#' @param y.column character. Name of the latitude column in `data` and `rdata`.
#' @param similarity.radii numeric vector of length two. Radii of influence in km: the first for the comparison with reference stations (Step 1) and the second for the comparison with validated nodes (Step 2). Default is `c(2, 0.2)`.
#' @param time.tolerance numeric or `NULL`. Time tolerance in minutes, for matching the timestamps of the tested node with the observations of the validated nodes in Step 2. If `NULL` (default), it is set to `3`. See \code{\link{find_neighbors}}.
#' @param similarity.duration numeric. Length of the rolling window, in hours. If one value is given, it is used for both steps; if two values are given, the first is used for Step 1 and the second for Step 2. Default is `c(336, 336)`.
#' @param compl.similarity.duration numeric. Minimum data completeness required within `similarity.duration`, expressed as a fraction of the hours in the window. If one value is given, it is used for both steps; if two values are given, the first is used for Step 1 and the second for Step 2. Default is `c(0.75, 0.75)`.
#' @param method character. Similarity indicator: `"nrmsd"` or `"mad"`. If one value is given, it is used for both steps; if two values are given, the first is used for Step 1 and the second for Step 2 (e.g., `c("nrmsd", "mad")`). Default is `"nrmsd"`. See \code{\link{spatial_similarity}}.
#' @param similarity.threshold numeric or `NULL`. Upper threshold of the similarity indicator. If one value is given, it is used for both steps; if two values are given, the first is used for Step 1 and the second for Step 2. If `NULL` (default), the threshold of each step is set according to its method: `0.68` for `"nrmsd"` and `5` for `"mad"`.
#' @param new logical. If `TRUE`, adds a column `<param>_chain_similarity` containing only the valid values; values with a flag other than `3` or `4` are replaced with `NA`. Default is `TRUE`.
#' @param param character or `NULL`. Base name of the new column (e.g., `"PM2.5"` gives `PM2.5_chain_similarity`). Used only if `new = TRUE`; if `NULL`, `data.column` is used. Default is `NULL`.
#'
#' @return The original dataset with:
#' \itemize{
#'   \item \code{flag_chain_similarity}: numeric. Quality flags based on the statistical test: \code{0} = insufficient data, \code{1} = isolated site, \code{2} = no similarity with the nearby higher-accuracy nodes,
#'   \code{3} = similar to a nearby validated node (indirect evaluation), \code{4} = similar to a nearby reference station (direct evaluation). It is `NA` where the tested value is missing.
#'   \item \code{<param>_chain_similarity}: numeric. Valid values only, with `NA` for the invalid ones (if `new = TRUE`).
#' }
#'
#' @export
chain_similarity_qctest <- function(data,
                                    date.column,
                                    data.column,
                                    id.column,
                                    rdata = NULL,
                                    rdata.column,
                                    rid.column,
                                    x.column,
                                    y.column,
                                    similarity.radii = c(2, 0.2),
                                    time.tolerance = NULL,
                                    similarity.duration = c(336, 336),
                                    compl.similarity.duration = c(0.75, 0.75),
                                    method = 'nrmsd',
                                    similarity.threshold = NULL,
                                    new = TRUE,
                                    param = NULL) {
  # new / param (checked before the early exit, which uses them)
  stopifnot(is.logical(new), length(new) == 1, !is.na(new))
  stopifnot(is.character(data.column), length(data.column) == 1, data.column %in% colnames(data))
  if (new) {
    if (is.null(param)) param <- data.column
    stopifnot(is.character(param), length(param) == 1, !is.na(param))
  }

  # Early exit if rdata is NULL
  if (is.null(rdata)) {
    data <- data %>% mutate(flag_chain_similarity = NA_real_)
    if (new) data <- data %>% mutate(!!paste0(param, "_chain_similarity") := NA_real_)
    return(data)
  }

  # arguments check
  stopifnot(is.character(date.column), date.column %in% colnames(data), date.column %in% colnames(rdata))
  stopifnot(is.character(id.column), length(id.column) == 1, id.column %in% colnames(data))
  stopifnot(is.character(rdata.column), rdata.column %in% colnames(rdata))
  stopifnot(is.character(rid.column), rid.column %in% colnames(rdata))
  stopifnot(is.character(x.column), x.column %in% colnames(data), x.column %in% colnames(rdata))
  stopifnot(is.character(y.column), y.column %in% colnames(data), y.column %in% colnames(rdata))

  # radii: exactly two values (Step 1, Step 2)
  stopifnot(is.numeric(similarity.radii), length(similarity.radii) == 2, all(similarity.radii > 0))

  # time tolerance in minutes (Step 2): NULL (default) or one non-negative value
  if (is.null(time.tolerance)) {
    time.tolerance <- 3
  }
  stopifnot(is.numeric(time.tolerance), length(time.tolerance) == 1, !is.na(time.tolerance), time.tolerance >= 0)

  # duration and completeness: one value (used for both steps) or two values
  stopifnot(is.numeric(similarity.duration), length(similarity.duration) %in% c(1, 2),
            all(similarity.duration > 0))
  stopifnot(is.numeric(compl.similarity.duration), length(compl.similarity.duration) %in% c(1, 2),
            all(compl.similarity.duration > 0), all(compl.similarity.duration <= 1))
  similarity.duration <- rep(similarity.duration, length.out = 2)
  compl.similarity.duration <- rep(compl.similarity.duration, length.out = 2)

  # method: one value (used for both steps) or two values
  stopifnot(is.character(method), length(method) %in% c(1, 2), all(method %in% c('mad', 'nrmsd')))
  method <- rep(method, length.out = 2)

  # threshold: NULL (default per method), one value or two values
  if (is.null(similarity.threshold)) {
    similarity.threshold <- ifelse(method == 'mad', 5, 0.68)
  } else {
    stopifnot(is.numeric(similarity.threshold), length(similarity.threshold) %in% c(1, 2),
              all(!is.na(similarity.threshold)), all(similarity.threshold > 0))
    similarity.threshold <- rep(similarity.threshold, length.out = 2)
  }
  #
  data <- data %>%
    mutate(!!id.column := as.character(.data[[id.column]])) %>%
    arrange(.data[[id.column]], .data[[date.column]], desc(!is.na(.data[[data.column]]))) %>%
    group_by(.data[[id.column]], .data[[date.column]]) %>%
    slice_head(n = 1) %>%
    ungroup()
  rdata <- rdata %>%
    mutate(!!rid.column := as.character(.data[[rid.column]])) %>%
    arrange(.data[[rid.column]], .data[[date.column]], desc(!is.na(.data[[rdata.column]]))) %>%
    group_by(.data[[rid.column]], .data[[date.column]]) %>%
    slice_head(n = 1) %>%
    ungroup()
  #
  # Part 1: Reference nodes
  reference_similarity_fun <- function(df) {
    # temporal alignment of the data and rdata
    rdata.agg <- rdata %>%
      group_by(across(all_of(c(rid.column, x.column, y.column)))) %>%
      arrange(.data[[date.column]]) %>%
      group_modify(~ {
        rdata.id <- .
        interp.values <- approx(x = rdata.id[[date.column]], y = rdata.id[[rdata.column]], xout = df[[date.column]], method = "linear")$y
        rdata.interp <- tibble(!!rid.column := rdata.id[[rid.column]][1], !!x.column := rdata.id[[x.column]][1], !!y.column := rdata.id[[y.column]][1],
                               !!date.column := df[[date.column]], !!rdata.column := interp.values)
        return(rdata.interp)}) %>%
      ungroup()
    #
    # find nearest neighbors (static references)
    neigh.info <- find_neighbors(data = df,
                                 date.column = date.column,
                                 id.column = id.column,
                                 rdata = rdata.agg,
                                 rdata.column = rdata.column,
                                 rid.column = rid.column,
                                 x.column = x.column,
                                 y.column = y.column,
                                 static = TRUE,
                                 match = 'exact',
                                 time.tolerance = 0,
                                 maxradii = similarity.radii[1],
                                 aggregation = FALSE)
    # no neighbors available
    if (is.null(neigh.info$data)) {
      df <- df %>% mutate(flag_chain_similarity = 1)
    } else {
      # calculate target - neighbor similarity
      target.data <- df %>%
        semi_join(neigh.info$neighbors %>% select(all_of(c(id.column, date.column))) %>% distinct(),
                  by = c(id.column, date.column)) %>%
        select(all_of(c(id.column, date.column, x.column, y.column, data.column))) %>%
        distinct()
      target.neigh.data <- neigh.info$data %>%
        select(all_of(c(rid.column, date.column, x.column, y.column, rdata.column)))
      target.neigh.similarity.results <- spatial_similarity(data = target.data,
                                                            date.column = date.column,
                                                            data.column = data.column,
                                                            id.column = id.column,
                                                            rdata = target.neigh.data,
                                                            rdata.column = rdata.column,
                                                            rid.column = rid.column,
                                                            similarity.duration = similarity.duration[1],
                                                            compl.similarity.duration = compl.similarity.duration[1],
                                                            similarity.threshold = similarity.threshold[1],
                                                            method = method[1]) %>%
        select(all_of(c(date.column, "flag_similarity"))) %>%
        rename(flag_chain_similarity = flag_similarity)
      df <- df %>%
        mutate(flag_chain_similarity = ifelse(is.na(.data[[data.column]]), NA, 1))
      df <- df %>%
        left_join(target.neigh.similarity.results, by = date.column, suffix = c("", "_new")) %>%
        mutate(flag_chain_similarity = coalesce(flag_chain_similarity_new, flag_chain_similarity)) %>%
        select(-flag_chain_similarity_new)
      df <- df %>%
        mutate(flag_chain_similarity = ifelse(is.na(.data[[data.column]]), NA_real_, flag_chain_similarity))
    }
    return(as_tibble(df))
  }
  data <- data %>% group_by(across(all_of(id.column))) %>%
    group_map(~ reference_similarity_fun(.x), .keep = TRUE) %>%
    bind_rows()
  # similarity with a reference station: 3 (spatial_similarity) -> 4 (direct evaluation)
  data <- data %>%
    mutate(flag_chain_similarity = ifelse(flag_chain_similarity == 3, 4, flag_chain_similarity))
  #
  # Part 2: higher accuracy nodes (not references)
  node_similarity_fun <- function(df) {
    id.target <- unique(df[[id.column]])
    df <- df %>%
      filter(flag_chain_similarity %in% c(0, 1))
    if (nrow(df) == 0) {
      df <- df %>%
        mutate(flag_chain_similarity_node = 0) %>%
        select(-flag_chain_similarity)
      return(df)
    }
    rdata <- data %>%
      filter(!(.data[[id.column]] %in% id.target), flag_chain_similarity == 4)
    if (nrow(rdata) == 0) {
      df <- df %>%
        mutate(flag_chain_similarity_node = 1) %>%
        select(-flag_chain_similarity)
      return(df)
    }
    #
    # find nearest neighbors (moving validated nodes, closest observation within the time tolerance)
    rdata.agg <- rdata %>%
      rename(!!paste0("r_", id.column) := all_of(id.column))
    neigh.info <- find_neighbors(data = df,
                                 date.column = date.column,
                                 id.column = id.column,
                                 rdata = rdata.agg,
                                 rid.column = paste0("r_", id.column),
                                 x.column = x.column,
                                 y.column = y.column,
                                 static = FALSE,
                                 match = 'near',
                                 time.tolerance = time.tolerance,
                                 maxradii = similarity.radii[2],
                                 aggregation = FALSE)
    # no neighbors available
    if (is.null(neigh.info$data)) {
      df <- df %>%
        mutate(flag_chain_similarity_node = 1) %>%
        select(-flag_chain_similarity)
      return(df)
    } else {
      # calculate target - neighbor similarity
      target.data <- df %>%
        filter(.data[[date.column]] %in% neigh.info$neighbors[[date.column]]) %>%
        select(all_of(c(id.column, date.column, x.column, y.column, data.column))) %>%
        distinct()
      # neighbor values at their matched observation, stamped with the target timestamp,
      # averaged over all matched neighbors into a single neighbor time series
      target.neigh.data <- neigh.info$neighbors %>%
        select(all_of(c(paste0("r_", id.column), date.column, paste0("r_", date.column)))) %>%
        distinct() %>%
        left_join(rdata.agg %>%
                    select(all_of(c(paste0("r_", id.column), date.column, data.column))) %>%
                    rename(!!paste0("r_", date.column) := all_of(date.column)),
                  by = c(paste0("r_", id.column), paste0("r_", date.column))) %>%
        group_by(across(all_of(date.column))) %>%
        summarise(!!data.column := mean(.data[[data.column]], na.rm = TRUE),
                  .groups = 'drop') %>%
        mutate(!!data.column := ifelse(is.nan(.data[[data.column]]), NA_real_, .data[[data.column]]),
               !!paste0("r_", id.column) := 'neighbors_mean') %>%
        select(all_of(c(paste0("r_", id.column), date.column, data.column)))
      target.neigh.similarity.results <- spatial_similarity(data = target.data,
                                                            date.column = date.column,
                                                            data.column = data.column,
                                                            id.column = id.column,
                                                            rdata = target.neigh.data,
                                                            rdata.column = data.column,
                                                            rid.column = paste0("r_", id.column),
                                                            similarity.duration = similarity.duration[2],
                                                            compl.similarity.duration = compl.similarity.duration[2],
                                                            similarity.threshold = similarity.threshold[2],
                                                            method = method[2]) %>%
        select(all_of(c(date.column, "flag_similarity"))) %>%
        rename(flag_chain_similarity = flag_similarity)
      #
      df <- df %>%
        left_join(target.neigh.similarity.results, by = date.column, suffix = c("_ref", "_node")) %>%
        mutate(flag_chain_similarity_node = ifelse(!is.na(flag_chain_similarity_node) & flag_chain_similarity_ref != 0, flag_chain_similarity_node, flag_chain_similarity_ref)) %>%
        select(-flag_chain_similarity_ref)
      df <- df %>%
        mutate(flag_chain_similarity_node = ifelse(is.na(.data[[data.column]]), NA_real_, flag_chain_similarity_node))
    }
    return(as_tibble(df))
  }
  data.nodes <- data %>%
    group_by(across(all_of(id.column))) %>%
    group_map(~ node_similarity_fun(.x), .keep = TRUE) %>%
    bind_rows()
  #
  data <- left_join(data, data.nodes %>% select(all_of(c(id.column, date.column, 'flag_chain_similarity_node'))), by = c(id.column, date.column))
  data <- data %>%
    mutate(flag_chain_similarity = case_when(flag_chain_similarity == 4 ~ 4,
                                             flag_chain_similarity == 2 ~ 2,
                                             flag_chain_similarity_node == 3 ~ 3,
                                             flag_chain_similarity %in% c(0, 1) & flag_chain_similarity_node == 2 ~ 2,
                                             flag_chain_similarity %in% c(0, 1) & (flag_chain_similarity_node == 1 | flag_chain_similarity == 1) ~ 1,
                                             flag_chain_similarity %in% c(0, 1) & (flag_chain_similarity_node == 0 | flag_chain_similarity == 0) ~ 0,
                                             TRUE ~ NA_real_)) %>%
    select(-flag_chain_similarity_node)
  #
  if (new) {
    data <- data %>%
      mutate(!!paste0(param, "_chain_similarity") := ifelse(flag_chain_similarity %in% c(3, 4), .data[[data.column]], NA_real_))
  }
  return(data)
}
