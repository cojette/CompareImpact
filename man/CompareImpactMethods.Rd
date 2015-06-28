\name{CompareImpactMethods}
\alias{CompareImpact}
\title{Creating a \code{CompareImpact} object}
\description{
  create the analysis results object returned by \code{CompareImpact()}.
}
\usage{
  CompareImpact(df, pre.period, post.period)
}
\author{
  JeongMin Kwon \email{cojette@gmail.com}
}
\examples{
\dontrun{
data(promo_data)
ci <- CompareImpact(promo_data, c(1,50), c(51,90))

}
}