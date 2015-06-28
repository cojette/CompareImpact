## CompareImpact Main Function
## Author: JeongMin Kwon <cojette@gmail.com>
## Note. Some inner functions are same as functions in CausalImpact package (https://github.com/google/CausalImpact)

devtools::install_github("google/CausalImpact")

.defaults <- list(niter = 1000,
                  standardize.data = TRUE,
                  prior.level.sd = 0.01,
                  nseasons = 1,
                  season.duration = 1,
                  dynamic.regression = FALSE)

# ------------------------------------------------------------------------------
FormatInputData <- function(data) {

  if (is.data.frame(data) && tolower(names(data)[1]) %in% c("date", "time")) {
    if (class(data$date) == "Date") {
      data <- zoo(data[, -1], data$date)
    } else {
      warning(paste0("Did you mean: data = zoo(data[, -1], data$",
                     names(data)[1], ")"))
    }
  }
  
  # Try to convert to zoo object
  data <- TryStop(as.zoo(data), "could not convert input data to zoo object")
  
  # Ensure <data> is formatted in such a way that rows represent time points
  if (is.null(ncol(data))) {
    dim(data) <- c(length(data), 1)
  }
  
  # Must have at least 3 time points
  assert_that(nrow(data) > 3)
  
  # Must not have NA in covariates (if any)
  if (ncol(data) >= 2) {
    assert_that(!any(is.na(data[, -1])))
  }
  return(data)
}

# ------------------------------------------------------------------------------
FormatInputPrePostPeriod <- function(pre.period, post.period, data) {
  # Checks and formats the <pre.period> and <post.period> input arguments.
  #
  # Args:
  # pre.period: two-element vector
  # post.period: two-element vector
  # data: already-checked zoo object, for reference only
  
  assert_that(!is.null(pre.period))
  assert_that(!is.null(post.period))
  assert_that(length(pre.period) == 2, length(post.period) == 2)
  assert_that(!any(is.na(pre.period)), !any(is.na(post.period)))
  if (class(time(data)) != class(pre.period) ||
        class(time(data)) != class(post.period)) {
    if (class(time(data)) == "integer") {
      pre.period <- as.integer(pre.period)
      post.period <- as.integer(post.period)
    } else if (class(time(data)) == "numeric") {
      pre.period <- as.numeric(pre.period)
      post.period <- as.numeric(post.period)
    } else {
      stop(paste0("pre.period (", class(pre.period), ") and post.period (",
                  class(post.period), ") should have the same class as the ",
                  "time points in the data (", class(time(data)), ")"))
    }
  }
  if (pre.period[1] < start(data)) {
    warning(paste0("Setting pre.period[1] to start of data: ", start(data)))
    pre.period[1] <- start(data)
  }
  if (pre.period[2] > end(data)) {
    warning(paste0("Setting pre.period[2] to end of data: ", end(data)))
    pre.period[2] <- end(data)
  }
  if (post.period[2] > end(data)) {
    warning(paste0("Setting post.period[2] to end of data: ", end(data)))
    post.period[2] <- end(data)
  }
  assert(pre.period[2] - pre.period[1] + 1 >= 3,
         "pre.period must span at least 3 time points")
  assert_that(post.period[2] >= post.period[1])
  assert_that(post.period[1] > pre.period[2])
  return(list(pre.period = pre.period, post.period = post.period))
}

# ------------------------------------------------------------------------------
FormatInputForCausalImpact <- function(data, pre.period, post.period,
                                       model.args, bsts.model,
                                       post.period.response, alpha) {
  # Checks and formats all input arguments supplied to CausalImpact(). See the
  # documentation of CausalImpact() for details.
  #

  assert(xor(!is.null(data) && !is.null(pre.period) && !is.null(post.period) &&
               is.null(bsts.model) && is.null(post.period.response),
             is.null(data) && is.null(pre.period) && is.null(post.period) &&
               !is.null(bsts.model) && !is.null(post.period.response)),
         paste0("must either provide data, pre.period, post.period, model.args",
                "; or bsts.model and post.period.response"))
  
  # Check <data> and convert to zoo, with rows representing time points
  if (!is.null(data)) {
    data <- FormatInputData(data)
  }
  
  # Check <pre.period> and <post.period>
  if (!is.null(data)) {
    checked <- FormatInputPrePostPeriod(pre.period, post.period, data)
    pre.period <- checked$pre.period
    post.period <- checked$post.period
  }
  
  # Parse <model.args>, fill gaps using <.defaults>
  model.args <- ParseArguments(model.args, .defaults)
  #
  # Check only those parts of <model.args> that are used in this file. The other
  # fields will be checked in FormatInputForConstructModel().
  #
  # Check <standardize.data>
  assert_that(length(model.args$standardize.data) == 1)
  assert_that(is.logical(model.args$standardize.data))
  assert_that(!is.na(model.args$standardize.data))
  
  # Check <bsts.model>
  if (!is.null(bsts.model)) {
    assert_that(class(bsts.model) == "bsts")
  }
  
  # Check <post.period.response>
  if (!is.null(bsts.model)) {
    assert_that(!is.null(post.period.response),
                is.vector(post.period.response),
                is.numeric(post.period.response))
  }
  
  # Check <alpha>
  assert_that(is.numeric(alpha))
  assert_that(length(alpha) == 1)
  assert_that(!is.na(alpha))
  assert_that(alpha > 0, alpha < 1)
  
  # Return updated arguments
  return(list(data = data, pre.period = pre.period, post.period = post.period,
              model.args = model.args, bsts.model = bsts.model,
              post.period.response = post.period.response, alpha = alpha))
}

# Main Function
#
# Args:
# data: Time series of response variable and any covariates. data.frame is recommended. 
# cbind(response vector, input vector1, input vector2....)
# Data from experiments with same environments is recommended
#
# pre.period: A vector of two indices specifying the first and the last
# time point of the pre-intervention period in the response
# vector \code{y}. 
#
# post.period: A vector of two indices specifying the first and the last day
# of the post-intervention period we wish to study. 
#
# model.args: Optional arguments that can be used to adjust the default
# construction of the state-space model used for inference.
# For full control over the model, you can construct your own
# model using the \code{bsts} package and feed the fitted model
# into \code{CausalImpact()} (see examples).
#
# bsts.model: Instead of passing in \code{data} and having
# \code{CausalImpact()} construct a model, it is possible to
# construct a model yourself using the \code{bsts} package. In
# this case, omit \code{data}, \code{pre.period}, and
# \code{post.period}. Instead only pass in \code{bsts.model},
# \code{y.post}, and \code{alpha} (optional). The model must
# have been fitted on data where the response variable was set
# to \code{NA} during the post-treatment period. The actual
# observed data during this period must then be passed to the
# function in \code{y.post}.
#
# post.period.response: Actual observed data during the post-intervention
# period. This is required if and only if a fitted
# \code{bsts.model} is passed instead of \code{data}.
#
# alpha: Desired tail-area probability for posterior intervals.
# Defaults to 0.05, which will produce central 95\% intervals.
#
# Returns:
# A CausalImpact object. This is a list of:
# series: observed data, counterfactual, pointwise and cumulative impact
# summary: summary table
# report: verbal description of the analysis
# model: contains bsts.model, the fitted model returned by bsts()
#
# Optional arguments for model.args:
# niter: number of MCMC iterations
# standardize.data: whether to standardize the data before model fitting
# prior.level.sd: standard deviation of the prior on the local level
# nseasons: number of seasons in the seasonal component
# season.duration: duration of each season
# dynamic.regression: whether to have dynamic instead of static coefficients
#
# return value: list of CausalImpact objects. You may use summary, plot functions in CausalImpact package individualy.

CompareImpact <- function(data = NULL,
                           pre.period = NULL,
                           post.period = NULL,
                           model.args = NULL,
                           bsts.model = NULL,
                           post.period.response = NULL,
                           alpha = 0.05) {

  
  data.length <- dim(data)[2]
  impact.list <- list()
  period.dim <- dim(pre.period)
  # Check input
  if (data.length <3 ) 
  {
    ## compare period
    if (!is.null(period.dim))
    {
      if(period.dim[2]!=2)
      {
        print(" Error: length(pre.period) not equal to 2 ")
        on.exit()
      }
      else 
      {
        for(j in 2:period.dim[1]) ## rbind data frame
        {
          impact.list[[j-1]] <- CausalImpact(data = data, 
                                             pre.period = pre.period[j-1,],
                                             post.period = pre.period[j,],
                                             model.args = model.args,
                                             bsts.model = bsts.model,
                                             post.period.response = post.period.response,
                                             alpha = alpha)
        }
      }
    }
    else 
    {
      print("Nothing to compare! You should use CausalImpact. ")
      on.exit()
    }
  }
  else
  {
    data.list <- list()
    for (i in 2:data.length)
    {
      data.list[[i]]<- cbind(data[,1], data[,i])
    }
    
    for (i in 2:data.length)
    {
      impact.list[[i-1]] <- CausalImpact(data = data.list[[i]], 
                                         pre.period = pre.period,
                                         post.period = post.period,
                                         model.args = model.args,
                                         bsts.model = bsts.model,
                                         post.period.response = post.period.response,
                                         alpha = alpha)
    }
  }
  
  ## return impact.list
  return(impact.list)
}

# ------------------------------------------------------------------------------
RunWithData <- function(data, pre.period, post.period, model.args, alpha) {

  # Zoom in on data in modeling range
  pre.period[1] <- max(pre.period[1], which(!is.na(data[, 1]))[1])
  data.modeling <- window(data, start = pre.period[1], end = post.period[2])
  if (is.null(ncol(data.modeling))) {
    dim(data.modeling) <- c(length(data.modeling), 1)
  }
  
  # Standardize all variables?
  UnStandardize <- identity
  if (model.args$standardize.data) {
    sd.results <- StandardizeAllVariables(data.modeling)
    data.modeling <- sd.results$data
    UnStandardize <- sd.results$UnStandardize
  }
  
  # Set observed response in post-period to NA
  window(data.modeling[, 1], start = post.period[1]) <- NA
  
  # Construct model and perform inference
  bsts.model <- ConstructModel(data.modeling, model.args)
  
  # Compile posterior inferences
  if (!is.null(bsts.model)) {
    y.post <- window(data[, 1], start = post.period[1], end = post.period[2])
    inferences <- CompilePosteriorInferences(bsts.model, y.post, alpha,
                                             UnStandardize)
  } else {
    inferences <- CompileNaInferences(data[, 1])
  }
  
  # Extend <series> to cover original range (padding with NA as necessary)
  empty <- zoo(, time(data))
  inferences$series <- merge(inferences$series, empty, all = TRUE)
  assert_that(nrow(inferences$series) == nrow(data))
  
  # Replace <y.model> by full original response
  inferences$series[, 1] <- data[, 1]
  
  # Assign response-variable names
  names(inferences$series)[1] <- "response"
  names(inferences$series)[2] <- "cum.response"
  
  # Return 'CausalImpact' object
  model <- list(pre.period = pre.period,
                post.period = post.period,
                model.args = model.args,
                bsts.model = bsts.model,
                alpha = alpha)
  impact <- list(series = inferences$series,
                 summary = inferences$summary,
                 report = inferences$report,
                 model = model)
  class(impact) <- "CausalImpact"
  return(impact)
}

CreateDataFrameForPlot <- function(impact) {
  # Creates a long-format data frame for CreateImpactPlot().
 
  # Create data frame from zoo series
  data <- as.data.frame(impact$series)
  data <- cbind(time = time(impact$series), data)
  
  # Reshape data frame
  tmp1 <- data[, c("time", "response", "point.pred", "point.pred.lower",
                   "point.pred.upper")]
  names(tmp1) <- c("time", "response", "mean", "lower", "upper")
  tmp1$baseline <- NA
  tmp1$metric <- "original"
  tmp2 <- data[, c("time", "response", "point.effect", "point.effect.lower",
                   "point.effect.upper")]
  names(tmp2) <- c("time", "response", "mean", "lower", "upper")
  tmp2$baseline <- 0
  tmp2$metric <- "pointwise"
  tmp2$response <- NA
  tmp3 <- data[, c("time", "response", "cum.effect", "cum.effect.lower",
                   "cum.effect.upper")]
  names(tmp3) <- c("time", "response", "mean", "lower", "upper")
  tmp3$metric <- "cumulative"
  tmp3$baseline <- 0
  tmp3$response <- NA
  data <- rbind(tmp1, tmp2, tmp3)
  data$metric <- factor(data$metric, c("original", "pointwise", "cumulative"))
  rownames(data) <- NULL
  return(data)
}

