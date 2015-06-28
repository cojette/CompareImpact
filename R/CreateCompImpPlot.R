## CompareImpact Plot Function
## Author: JeongMin Kwon <cojette@gmail.com>
## Note. Some inner functions are same as functions in CausalImpact package (https://github.com/google/CausalImpact)

CreatePeriodMarkers <- function(pre.period, post.period, time.range) {
  # Creates a vector of period markers to display.
  #
  # Args:
  #   pre.period: vector of 2 time points that define the pre-period
  #   post.period: vector of 2 time points that define the post-period
  #   time.range: vector of 2 elements specifying range of timepoints in data
  #
  # Returns:
  #   vector of period markers that should be displayed
  
  idx <- NULL
  if (pre.period[1] > time.range[1]) {
    idx <- c(idx, pre.period[1])
  }
  if (pre.period[2] < post.period[1] - 1) {
    idx <- c(idx, pre.period[2])
  }
  idx <- c(idx, post.period[1] - 1)
  if (post.period[2] < time.range[2]) {
    idx <- c(idx, post.period[2])
  }
  class(idx) <- class(pre.period)
  return(as.numeric(idx))
}


CreateCompImpPlot <- function(impact.list, metrics = c("original", "pointwise",
                                                       "cumulative")) {
  # Creates a plot of observed data and counterfactual predictions.
  #
  # Args:
  #   impact.list:  list of \code{CausalImpact} results object returned by
  #            \code{CausalImpact()}.
  #   metrics: Which metrics to include in the plot. Can be any combination of
  #            "original", "pointwise", and "cumulative". Default value is all. 
  #
  # Returns:
  #   A ggplot2 object that can be plotted 
  imp.length <- length(impact.list)
  data <- list()
  
  metrics <- sapply(metrics, function(m) match.arg(m, c("original", "pointwise",
                                                        "cumulative")))
  q <- plot(impact.list[[1]])
  
  ## first case: various timeline
  if (identical(impact.list[[2]]$model$pre.period,impact.list[[1]]$model$pre.period))
  {
    #     q <- plot(impact.list[[1]])
    
    for(i in 2:length(impact.list))
    {
      if(impact.list[[i]]$summary$p[1]>impact.list[[i]]$summary$alpha[1])
      {
        print(paste0("series ", i, " would generally not be considered statistically significant and plot is omitted."))
      }
      else
      {
        data <- CreateDataFrameForPlot(impact.list[[i]])
        metrics <- sapply(metrics, function(m) match.arg(m, c("original", "pointwise",
                                                              "cumulative")))
        data <- data[data$metric %in% metrics, ]
        data$metric <- factor(data$metric, metrics)      
        
        q <- q + geom_ribbon(aes(ymin = lower, ymax = upper),
                             as.vector(data), fill = i, alpha = 0.2)
        q <- q + geom_line(aes(y = mean), data,
                           size = 0.6, colour = i, linetype = "dashed")
      }
    }
  }
  
  ## Second: various intercept
  else
  {
    for(i in 2:length(impact.list))
    {
      if(impact.list[[i]]$summary$p[1]>impact.list[[i]]$summary$alpha[1])
      {
        print(paste0("series ", i, " would generally not be considered statistically significant and plot is omitted."))
      }
      else
      {
        data <- CreateDataFrameForPlot(impact.list[[i]])
        metrics <- sapply(metrics, function(m) match.arg(m, c("original", "pointwise",
                                                              "cumulative")))
        data <- data[data$metric %in% metrics, ]
        data$metric <- factor(data$metric, metrics)      
        
        xintercept <- CreatePeriodMarkers(impact.list[[1]]$model$pre.period,
                                          impact.list[[i]]$model$post.period,
                                          range(data$t))
        q <- q + geom_vline(xintercept = xintercept,
                            colour = "darkgrey", size = 0.8, linetype = "dashed")
        
        q <- q + geom_ribbon(aes(ymin = lower, ymax = upper),
                             as.vector(data), fill = i, alpha = 0.2)
        q <- q + geom_line(aes(y = mean), data,
                           size = 0.6, colour = i, linetype = "dashed")
      }
    }
  }

  return(q)
}


