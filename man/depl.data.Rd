\name{depl.data}
\alias{depl.data}
\docType{data}
\title{
Depletion data
%%   ~~ data name/kind ... ~~
}
\description{
Counts of fish observed from a mulitple pass removal method.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(depl.data)}
\format{
  A data frame with 166 observations on the following 8 variables.
  \describe{
    \item{\code{site.name}}{a factor with levels \code{Big0Springs-1} \code{Big0Springs-4} \code{Big0Springs-5} \code{Big0Springs-8} \code{BRM-223986} \code{CMP-030730} \code{CMP-330226} \code{CMP-383986} \code{DFC-438922} \code{DUN-290498} \code{GBD-189938} \code{LEM-000002} \code{LEM-000003} \code{LEM-00001B} \code{LEM-00001C} \code{LEM-00001D} \code{LEM-000943} \code{LEM-001487} \code{LEM-002095} \code{LEM-004959} \code{LEM-005583} \code{LEM-006575} \code{LEM-007007} \code{LEM-009055} \code{LEM-009135} \code{LEM-010703} \code{LEM-013151} \code{LEM-014799} \code{LEM-019295} \code{LEM-025519} \code{LEM-026031} \code{LEM-027487} \code{LEM-028079} \code{LEM-029103} \code{LEM-029535} \code{LEM-030559} \code{LEM-034767} \code{LEM-035679} \code{LEM-036015} \code{LEM-038111} \code{LEM-038351} \code{LEM-039375} \code{LEM-040911} \code{LEM-042447} \code{LEM-043471} \code{LEM-045487} \code{LEM-047023} \code{LEM-047567} \code{LEM-048559} \code{LEM-048847} \code{LEM-049615} \code{LEM-053087} \code{LEM-058207} \code{LEM-060847} \code{LEM-061903} \code{LEM-063951} \code{LEM-064351} \code{LEM-067023} \code{LEM-067535} \code{LEM-068047} \code{LEM-076751} \code{LEM-078287} \code{LEM-080335} \code{LEM-083807} \code{LEM-084831} \code{LEM-086063} \code{LEM-087503} \code{LEM-091599} \code{LEM-093023} \code{LEM-100191} \code{LEM-100303} \code{LEM-101551} \code{LEM-102239} \code{LEM-107439} \code{LEM-109999} \code{LEM-111455} \code{LEM-119215} \code{LEM-121695} \code{LEM-128175} \code{LEM-130479} \code{LEM-136031} \code{LEM-141151} \code{LEM-149967} \code{LEM-156111} \code{LEM-170447} \code{LEM-193375} \code{LEM-197071} \code{LEM-207711} \code{LEM-217551} \code{LEM-222671} \code{LEM-240479} \code{LEM-252767} \code{LEM-258911} \code{LEM-275935} \code{LEM-306015} \code{LEM-315823} \code{LEM-358863} \code{LEM-362335} \code{LEM-394703} \code{LEM-395727} \code{LEM-413535} \code{LEM-423775} \code{LEM-442831} \code{LEM-449375} \code{LEM-452047} \code{LEM-455519} \code{LEM-482143} \code{LEM-510383} \code{Little0Springs2} \code{MFJD-207602} \code{MFJD-275954} \code{MFJD-289522} \code{MFJD-290034} \code{MFJD-298738} \code{MFJD-484850} \code{MFJD-497650} \code{MRC-011122} \code{MRC-018034} \code{MRC-067186} \code{MRC-095602} \code{MRC-183154} \code{MRC-214642} \code{MRC-232818} \code{SFM-292210} \code{SFM-396146} \code{SFS-000287} \code{SFS-010642} \code{SFS-025887} \code{SFS-030354} \code{SFS-045010} \code{SFS-055583} \code{SFS-073682} \code{SFS-079506} \code{SFS-087698} \code{SFS-096543} \code{SFS-115999} \code{SFS-159007} \code{SFS-248210} \code{SFS-307154} \code{SFS-339922} \code{SFS-348959} \code{SFS-382610} \code{SFS-404767} \code{SFS-428434} \code{SUM-003826} \code{SUM-325362} \code{SUM-429810} \code{TEX-000066} \code{TEX-167282} \code{TEX-429426} \code{VGR-144114} \code{VGR-404210}}
    \item{\code{year}}{a character vector}
    \item{\code{stream.name}}{a factor with levels \code{Agency Creek} \code{Basin Creek} \code{Bear Valley Creek} \code{Big Bear Creek} \code{Big Eightmile Creek} \code{Big Springs Creek} \code{Big Timber Creek} \code{Bohannon Creek} \code{Bridge Creek} \code{Buck Creek} \code{Camp Creek} \code{Canyon Creek} \code{Cruikshank Creek} \code{Dry Fork Clear Creek} \code{Duncan Creek} \code{East Fork Bohannon Creek} \code{East Fork Hayden Creek} \code{East Fork South Fork Salmon River} \code{Eighteenmile Creek} \code{Flume Creek} \code{Granite Boulder Creek} \code{Grouse Creek} \code{Hawley Creek} \code{Hayden Creek} \code{Kenney Creek} \code{Lake Creek} \code{Lemhi Little Springs Creek} \code{Lemhi River} \code{Lick Creek} \code{Little Eightmile Creek} \code{Little Springs} \code{Little Springs Creek} \code{Little Timber Creek} \code{McDevitt Creek} \code{Middle Fork John Day River} \code{Murderers Creek} \code{Pattee Creek} \code{Reservoir Creek} \code{Secesh River} \code{South Fork Murderers Creek} \code{Summit Creek} \code{Tex Creek} \code{Texas Creek} \code{Trail Creek} \code{Unnamed Tributary to Lower Lemhi River} \code{Unnamed Tributary to Upper Lemhi River} \code{Vinegar Creek} \code{West Fork Wimpey Creek} \code{Wimpey Creek} \code{Withington Creek} \code{Wright Creek} \code{Zena Creek}}
    \item{\code{n.pass}}{a numeric vector}
    \item{\code{pass.1}}{a numeric vector}
    \item{\code{pass.2}}{a numeric vector}
    \item{\code{pass.3}}{a numeric vector}
    \item{\code{pass.4}}{a numeric vector}
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
data(depl.data)
## maybe str(depl.data) ; plot(depl.data) ...
}
\keyword{datasets}
