#' Quality Control based on mFILTER
#'
#' Applies the full quality control (QC) procedure of mFILTER in the raw data of mobile sensor.
#'
#' @details
#' This function performs all quality control (QC) algorithms included in the mFILTER approach for mobile sensors. QC is applied using a six-step procedure in the following order: \cr\cr
#'
#' \strong{Raw Data} ---------  \strong{Location} --> \strong{Range Validity} --> \strong{Constant Value} --> \strong{Spatial Correlation} --> \strong{Spatial Similarity} --> \strong{Chain Similarity} ---------  \strong{QC Raw Data}\cr\cr
#'
#' The \code{\link{location_qctest}}, \code{\link{range_validity_qctest}}, and \code{\link{constant_value_qctest}} tests are applied directly on the raw data.
#'
#' The \code{\link{spatial_correlation_qctest}} is applied only on the measurements that passed the first three tests
#' (\code{flag_location = 1}, \code{flag_range = 1} and \code{flag_constant_value = 2}).
#'
#' The \code{\link{spatial_similarity_qctest}} and \code{\link{chain_similarity_qctest}} are applied only on the measurements that additionally passed the
#' \code{\link{spatial_correlation_qctest}}, namely \code{flag_location = 1}, \code{flag_range = 1}, \code{flag_constant_value = 2} and \code{flag_correlation = 3}.
#'
#' The measurements that did not pass a previous step are set to \code{NA} before the next step, so that they are
#' neither tested nor used as neighbors. At the end of the process, the flags of all tests are included in the dataset.
#'
#' @param data A data frame containing the dataset to be tested.
#' @param date.column character. Name of the date/time column in `data`.
#' @param data.column character. Name of the column in `data` on which the test is performed.
#' @param id.column character. Name or ID of the measurement location (node) in `data`.
#' @param x.column see in \code{\link{location_qctest}}, \code{\link{spatial_correlation_qctest}}, \code{\link{spatial_similarity_qctest}}, and \code{\link{chain_similarity_qctest}}.
#' @param y.column see in \code{\link{location_qctest}}, \code{\link{spatial_correlation_qctest}}, \code{\link{spatial_similarity_qctest}}, and \code{\link{chain_similarity_qctest}}.
#' @param v.th see in \code{\link{location_qctest}}.
#' @param upper.bound see in \code{\link{range_validity_qctest}}.
#' @param lower.bound see in \code{\link{range_validity_qctest}}.
#' @param persist.duration see in \code{\link{constant_value_qctest}}.
#' @param compl.persist.duration see in \code{\link{constant_value_qctest}}.
#' @param min.variation see in \code{\link{constant_value_qctest}}.
#' @param rdata see in \code{\link{spatial_correlation_qctest}}, \code{\link{spatial_similarity_qctest}}, and \code{\link{chain_similarity_qctest}}.
#' @param rdata.column see in \code{\link{spatial_correlation_qctest}}, \code{\link{spatial_similarity_qctest}}, and \code{\link{chain_similarity_qctest}}.
#' @param rid.column see in \code{\link{spatial_correlation_qctest}}, \code{\link{spatial_similarity_qctest}}, and \code{\link{chain_similarity_qctest}}.
#' @param compl.corr.duration see in \code{\link{spatial_correlation_qctest}}.
#' @param corr.radii see in \code{\link{spatial_correlation_qctest}}.
#' @param corr.method see in \code{\link{spatial_correlation_qctest}}.
#' @param corr.th see in \code{\link{spatial_correlation_qctest}}.
#' @param historical see in \code{\link{spatial_correlation_qctest}}.
#' @param ndays see in \code{\link{spatial_correlation_qctest}}.
#' @param similarity.radii `similarity.radii` in \code{\link{spatial_similarity_qctest}}.
#' @param similarity.duration `similarity.duration` in \code{\link{spatial_similarity_qctest}}.
#' @param compl.similarity.duration `compl.similarity.duration` in \code{\link{spatial_similarity_qctest}}.
#' @param similarity.method `method` in \code{\link{spatial_similarity_qctest}}.
#' @param similarity.threshold `similarity.threshold` in \code{\link{spatial_similarity_qctest}}.
#' @param chain.radii `similarity.radii` in \code{\link{chain_similarity_qctest}}.
#' @param chain.duration `similarity.duration` in \code{\link{chain_similarity_qctest}}.
#' @param compl.chain.duration `compl.similarity.duration` in \code{\link{chain_similarity_qctest}}.
#' @param chain.method `method` in \code{\link{chain_similarity_qctest}}.
#' @param chain.threshold `similarity.threshold` in \code{\link{chain_similarity_qctest}}.
#' @param param see in \code{\link{range_validity_qctest}}, \code{\link{constant_value_qctest}}, \code{\link{spatial_correlation_qctest}}, \code{\link{spatial_similarity_qctest}}, and \code{\link{chain_similarity_qctest}}.
#' @param new see in \code{\link{range_validity_qctest}}, \code{\link{constant_value_qctest}}, \code{\link{spatial_correlation_qctest}}, \code{\link{spatial_similarity_qctest}}, and \code{\link{chain_similarity_qctest}}.
#'
#' @return Returns the original data including the quality flags of all tests for each measurement.
#'
#' @seealso \code{\link{location_qctest}}, \code{\link{range_validity_qctest}}, \code{\link{constant_value_qctest}}, \code{\link{spatial_correlation_qctest}}, \code{\link{spatial_similarity_qctest}}, and \code{\link{chain_similarity_qctest}}
#'
#' @export
run_mFILTER_QC <- function(data,
                           date.column,
                           data.column,
                           id.column,
                           x.column,
                           y.column,
                           rdata = NULL,
                           rdata.column = NULL,
                           rid.column = NULL,
                           # location test arguments
                           v.th = 140,
                           # range test arguments
                           upper.bound = 1000,
                           lower.bound = 0,
                           # constant value test arguments
                           persist.duration = 3,
                           compl.persist.duration = 0.6,
                           min.variation = 0.1,
                           # spatial correlation arguments
                           compl.corr.duration = 0.6,
                           corr.radii = 10,
                           corr.method = 'lut',
                           corr.th = 0.75,
                           historical = TRUE,
                           ndays = 30,
                           # spatial similarity arguments
                           similarity.radii = 30,
                           similarity.duration = 72,
                           compl.similarity.duration = 0.6,
                           similarity.method = 'mad',
                           similarity.threshold = 5,
                           # chain similarity arguments
                           chain.radii = c(2, 0.2),
                           chain.duration = c(336, 336),
                           compl.chain.duration = c(0.1, 0.1),
                           chain.method = 'nrmsd',
                           chain.threshold = 0.68,
                           # other arguments
                           param = 'PM2.5',
                           new = FALSE) {
  #
  # temporary columns with the measurements that passed the previous steps
  flagged.column <- paste0(data.column, '_flagged')
  correlation.column <- paste0(data.column, '_correlation')
  #
  # mFILTER approach (Mobile Sensors)
  data <- location_qctest(data, date.column = date.column, id.column = id.column, x.column = x.column, y.column = y.column, v.th = v.th)
  data <- range_validity_qctest(data, data.column = data.column, upper.bound = upper.bound, lower.bound = lower.bound, param = param, new = new)
  data <- constant_value_qctest(data, date.column = date.column, data.column = data.column, id.column = id.column, persist.duration = persist.duration,
                                compl.persist.duration = compl.persist.duration, min.variation = min.variation, param = param, new = new)
  # measurements valid in location, range and constant value tests
  data <- data %>%
    mutate(!!paste0(data.column, '_flagged') := ifelse(.data[['flag_location']] == 1 & .data[['flag_range']] == 1 & .data[['flag_constant_value']] == 2,
                                      .data[[data.column]], NA))
  data <- spatial_correlation_qctest(data, date.column = date.column, data.column = paste0(data.column, '_flagged'), id.column = id.column,
                                     rdata = rdata, rdata.column = rdata.column, rid.column = rid.column, x.column = x.column, y.column = y.column,
                                     compl.corr.duration = compl.corr.duration, corr.radii = corr.radii, corr.method = corr.method, corr.th = corr.th,
                                     historical = historical, ndays = ndays, param = param, new = new)
  # measurements valid also in the spatial correlation test
  data <- data %>%
    mutate(!!paste0(data.column, '_correlation') := ifelse(!is.na(.data[[paste0(data.column, '_flagged')]]) & .data[['flag_correlation']] == 3,
                                          .data[[data.column]], NA))
  data <- spatial_similarity_qctest(data, date.column = date.column, data.column = paste0(data.column, '_correlation'), id.column = id.column,
                                    rdata = rdata, rdata.column = rdata.column, rid.column = rid.column, x.column = x.column, y.column = y.column,
                                    similarity.radii = similarity.radii, similarity.duration = similarity.duration,
                                    compl.similarity.duration = compl.similarity.duration, method = similarity.method,
                                    similarity.threshold = similarity.threshold, param = param, new = new)
  data <- chain_similarity_qctest(data, date.column = date.column, data.column = paste0(data.column, '_correlation'), id.column = id.column,
                                  rdata = rdata, rdata.column = rdata.column, rid.column = rid.column, x.column = x.column, y.column = y.column,
                                  similarity.radii = chain.radii, similarity.duration = chain.duration,
                                  compl.similarity.duration = compl.chain.duration, method = chain.method,
                                  similarity.threshold = chain.threshold, param = param, new = new)
  # removal of the temporary columns
  data <- data %>%
    select(-all_of(c(paste0(data.column, '_flagged'), paste0(data.column, '_correlation'))))
  #
  return(data)
}
