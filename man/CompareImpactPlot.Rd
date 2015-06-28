\name{CompareImpactPlot}
\alias{CreateCompImpPlot}
\title{Plotting a \code{CompareImpact} object}
\description{
  Methods for plotting the results of an
  analysis results object returned by \code{CompareImpact()}.
}
\usage{
  CreateCompImpPlot(object,...)
}
\arguments{
  \item{object}{A \code{CompareImpact()} results object, as returned
  by \code{CompareImpact()}.}

   \item{...}{Optional additional arguments. You can choose one option original, pointwise, cumulative, and by default, all three options are plotted. }
}
\author{
  JeongMin Kwon \email{cojette@gmail.com}
}
\examples{
\dontrun{
ci <- CompareImpact(...)

# Create a plot
CreateCompImpPlot(ci)

# Customize a plot
ci.plot <- CreateCompImpPlot(ci)
ci.plot <- ci.plot + theme_bw(base_size = 20)
}
}