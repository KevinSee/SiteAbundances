\name{MR.data}
\alias{MR.data}
\docType{data}
\title{
Mark-Recapture data
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
Counts of fish from a 2-pass mark recapture methodology. M is the number of fish caught on the first pass, all of whom are marked. C is the number of fish caught on the second pass. R is the number of fish caught on both passes. 
}
\usage{data(MR.data)}
\format{
  A data frame with 93 observations on the following 5 variables.
  \describe{
    \item{\code{SiteName}}{a factor with levels \code{Beyeler-1} \code{Beyeler-9} \code{Big0Springs-5} \code{Big0Springs-6} \code{CanyonMega1} \code{CanyonMega2} \code{CanyonMega4} \code{CanyonMega5} \code{CanyonMega6} \code{Emmigration} \code{KennyMega1} \code{KennyMega2} \code{KennyMega3} \code{LEM-001487} \code{LEM-006575} \code{LEM-020943} \code{LEM-029103} \code{LEM-029135} \code{LEM-029535} \code{LEM-038111} \code{LEM-045487} \code{LEM-045519} \code{LEM-048559} \code{LEM-048847} \code{LEM-053711} \code{LEM-054495} \code{LEM-061903} \code{LEM-091567} \code{LEM-119215} \code{LEM-160207} \code{LEM-190895} \code{LEM-223663} \code{LemhiMega1} \code{LemhiMega2} \code{LemhiMega3} \code{LemhiMega4} \code{LemhiMega5} \code{LemhiMega6} \code{Lick Cr-1} \code{Lick Cr-2} \code{Secesh Mega 1} \code{Secesh Mega 2} \code{Secesh Mega 3} \code{Secesh Mega 4} \code{Secesh Mega 5} \code{Secesh Mega 6} \code{SFS-000287} \code{SFS-017695} \code{SFS-025887} \code{SFS-029471} \code{SFS-073682} \code{SFS-079506} \code{SFS-086815} \code{SFS-115999} \code{SFS-159007} \code{SFS-176082} \code{SFS-181535} \code{SFS-206111} \code{SFS-217887} \code{SFS-231826} \code{SFS-235154} \code{SFS-241618} \code{SFS-279839} \code{SFS-308882} \code{SFS-312607} \code{SFS-317074} \code{SFS-339922} \code{SFS-340767} \code{SFS-345375} \code{SFS-348959} \code{SFS-374418} \code{SFS-382610} \code{SFS-404767} \code{SFS-412050} \code{SFS-428434} \code{SFS-444818} \code{SFS-488223} \code{SFS-505490}}
    \item{\code{Year}}{a factor with levels \code{2009} \code{2010} \code{2011} \code{2012}}
    \item{\code{M}}{a numeric vector}
    \item{\code{C}}{a numeric vector}
    \item{\code{R}}{a numeric vector}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(MR.data)
## maybe str(MR.data) ; plot(MR.data) ...
}
\keyword{datasets}
