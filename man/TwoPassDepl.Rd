\name{TwoPassDepl}
\alias{TwoPassDepl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Two Pass Depletion Abundance Estimator
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
a function to estimate abundance based on a two pass depletion method
}
\usage{
TwoPassDepl(data, rmInvalid=FALSE, p.threshold=0.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ dataframe with two columns containing the counts from each pass }
  \item{rmInvalid}{ if TRUE, make any abundances NA where the probability of capture is estimated to be below p.threshold}
  \item{p.threshold}{ the user-defined cutoff for declaring an estimate invalid. Set to 0.2 as a default}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
    \item{N.hat}{ - estimate of abundance}
    \item{N.hat.SE}{ - standard error of the estimate of abundance}
    \item{p.hat}{ - estimate of the probability of capture}
}
\references{
Lockwood, Roger N. and J. C. Schneider. 2000. Stream fish population estimates by markand-recapture and depletion methods. Chapter 7 in Schneider, James C. (ed.) 2000. Manual of fisheries survey methods II: with periodic updates. Michigan Department of Natural Resources, Fisheries Special Report 25, Ann Arbor.
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
# pull out examples with two-pass depletion data
two.pass.data = subset(depl.data, n.pass==2, select=c('pass.1', 'pass.2'))
# estimate abundance
TwoPassDepl(two.pass.data)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
