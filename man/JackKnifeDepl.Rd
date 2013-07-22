\name{JackKnifeDepl}
\alias{JackKnifeDepl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Jacknife Removal Method Estimator
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A function to estimate abundance based on a removal (depletion) method. Particularly useful when the number of total passes is small (<5), the total catch is small (<30) and the capture probability is low (<0.4).
}
\usage{
JackKnifeDepl(data, group = NULL, grouped = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{dataframe or matrix with three columns containing the counts from each pass.}
  \item{group}{factor describing which group that depletion comes from. If grouped==FALSE, not needed}
  \item{grouped}{logical. If FALSE, all sites are treated as having the same probability of capture. If TRUE, each group has a probability of capture estimated separately.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{N.hat}{ - estimate of abundance}
  \item{p.hat}{ - estimate of the probability of capture}
}
\references{
Hankin, D. and Mohr, M. (2001). Improved two-phase survey designs for estimation of fish abundance in small streams. Unpublished manuscript.
}
\author{
Kevin See, QCI, Seattle, WA
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(depl.data)
# pull out examples with three-pass depletion data
three.pass.data = subset(depl.data, n.pass==3, select=c('pass.1', 'pass.2', 'pass.3', 'stream.name'))
# estimate abundance, ungrouped data
JackKnifeDepl(three.pass.data[,1:3])

# estimate abundance, group data by stream name
JackKnifeDepl(three.pass.data[,1:3], group=three.pass.data$stream.name, grouped=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
