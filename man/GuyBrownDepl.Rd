\name{GuyBrownDepl}
\alias{GuyBrownDepl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
3 pass depletion estimator
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
a function to estimate abundance based on a three pass depletion method
}
\usage{
GuyBrownDepl(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{dataframe or matrix with three columns containing the counts from each pass.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a dataframe consisting of 3 columns:
  \item{N.hat}{ - estimate of abundance}
  \item{N.hat.SE}{ - standard error of the estimate of abundance}
  \item{p.hat}{ - estimate of the probability of capture}

}
\references{
Hayes, D.B., J.R. Bence, T.J. Kwak, and B.E. Thompson. 2007. Abundance, biomass, and production. Pages 327-374 in C.S. Guy and M.L. Brown, editors. Analysis and interpretation of freshwater fisheries data. American Fisheries Society, Bethesda, Maryland.
%% ~put references to the literature/web site here ~
}
\author{
Kevin See, QCI, Seattle, WA
%%  ~~who you are~~
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
three.pass.data = subset(depl.data, n.pass==3, select=c('pass.1', 'pass.2', 'pass.3'))
# estimate abundance
GuyBrownDepl(three.pass.data)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
