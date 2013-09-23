\name{fish.data}
\alias{fish.data}
\docType{data}
\title{
ISEMP fish data
}
\description{
Fish data from mark recapture, three pass depletion and single pass surveys of streams.
}
\usage{data(fish.data)}
\format{
  A data frame with 627 observations on the following 8 variables.
  \describe{
    \item{\code{Site.ID}}{a factor with levels \code{1} \code{2} \code{3} \code{4} \code{5} \code{6} \code{7} \code{8} \code{9} \code{10} \code{11} \code{12} \code{13} \code{14} \code{15} \code{16} \code{17} \code{18} \code{19} \code{20} \code{21} \code{22} \code{23} \code{24} \code{25} \code{26} \code{27} \code{28} \code{29} \code{30} \code{31} \code{32} \code{33} \code{34} \code{35} \code{36} \code{37} \code{38} \code{39} \code{40} \code{41} \code{42} \code{43} \code{44} \code{45} \code{46} \code{47} \code{48} \code{49} \code{50} \code{51} \code{52} \code{53} \code{54} \code{55} \code{56} \code{57} \code{58} \code{59} \code{60} \code{61} \code{62} \code{63} \code{64} \code{65} \code{66} \code{67} \code{68} \code{69} \code{70} \code{71} \code{72} \code{73} \code{74} \code{75} \code{76} \code{77} \code{78} \code{79} \code{80} \code{81} \code{82} \code{83} \code{84} \code{85} \code{86} \code{87} \code{88} \code{89} \code{90} \code{91} \code{92} \code{93} \code{94} \code{95} \code{96} \code{97} \code{98} \code{99} \code{100} \code{101} \code{102} \code{103} \code{104} \code{105} \code{106} \code{107} \code{108} \code{109} \code{110} \code{111} \code{112} \code{113} \code{114} \code{115} \code{116} \code{117} \code{118} \code{119} \code{120} \code{121} \code{122} \code{123} \code{124} \code{125} \code{126} \code{127} \code{128} \code{129} \code{130} \code{131} \code{132} \code{133} \code{134} \code{135} \code{136} \code{137} \code{138} \code{139} \code{140} \code{141} \code{142} \code{143} \code{144} \code{145} \code{146} \code{147} \code{148} \code{149} \code{150} \code{151} \code{152} \code{153} \code{154} \code{155} \code{156} \code{157} \code{158} \code{159} \code{160} \code{161} \code{162} \code{163} \code{164} \code{165} \code{166} \code{167} \code{168} \code{169} \code{170} \code{171} \code{172} \code{173} \code{174} \code{175} \code{176} \code{177} \code{178} \code{179} \code{180} \code{181} \code{182} \code{183} \code{184} \code{185} \code{186} \code{187} \code{188} \code{189} \code{190} \code{191} \code{192} \code{193} \code{194} \code{195} \code{196} \code{197} \code{198} \code{199} \code{200} \code{201} \code{202} \code{203} \code{204} \code{205} \code{206} \code{207} \code{208} \code{209} \code{210} \code{211} \code{212} \code{213} \code{214} \code{215} \code{216} \code{217} \code{218} \code{219} \code{220} \code{221} \code{222} \code{223} \code{224} \code{225} \code{226} \code{227} \code{228} \code{229} \code{230} \code{231} \code{232} \code{233} \code{234} \code{235} \code{236} \code{237} \code{238} \code{239} \code{240} \code{241} \code{242} \code{243} \code{244} \code{245} \code{246} \code{247} \code{248} \code{249} \code{250} \code{251} \code{252} \code{253} \code{254} \code{255} \code{256} \code{257} \code{258} \code{259} \code{260} \code{261} \code{262} \code{263} \code{264} \code{265} \code{266} \code{267} \code{268} \code{269} \code{270} \code{271} \code{272} \code{273} \code{274} \code{275} \code{276} \code{277} \code{278} \code{279} \code{280} \code{281} \code{282} \code{283} \code{284} \code{285} \code{286} \code{287} \code{288} \code{289} \code{290} \code{291} \code{292} \code{293} \code{294} \code{295} \code{296} \code{297} \code{298} \code{299} \code{300} \code{301} \code{302} \code{303} \code{304} \code{305} \code{306} \code{307} \code{308} \code{309} \code{310} \code{311} \code{312} \code{313} \code{314} \code{315} \code{316} \code{317} \code{318} \code{319} \code{320} \code{321} \code{322} \code{323} \code{324} \code{325} \code{326} \code{327} \code{328} \code{329} \code{330} \code{331} \code{332} \code{333} \code{334} \code{335} \code{336} \code{337} \code{338} \code{339} \code{340} \code{341} \code{342} \code{343} \code{344} \code{345} \code{346}}
    \item{\code{Subbasin}}{a factor with levels \code{A} \code{B} \code{C} \code{D} \code{E} \code{F} \code{G} \code{H}}
    \item{\code{Species}}{a factor with levels \code{Chinook Salmon} \code{Steelhead Trout}}
    \item{\code{Method}}{a factor with levels \code{Depletion} \code{Mark Recapture} \code{Single Pass} \code{Snorkel}}
    \item{\code{SiteLength}}{a numeric vector}
    \item{\code{Pass1.M}}{a numeric vector}
    \item{\code{Pass2.C}}{a numeric vector}
    \item{\code{Pass3.R}}{a numeric vector}
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
data(fish.data)
## maybe str(fish.data) ; plot(fish.data) ...
}
\keyword{datasets}
