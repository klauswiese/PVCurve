\name{PVCurveMean}
\alias{PVCurveMean}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Preassurve - Volume Curve with estimete of Osmosis Cero Point by the mean of the two points where the rupture is observe
}
\description{
Estimate parameteres based on Preassurve - Volume Curve}
\usage{
PVCurve(x,y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{1 - Contenido Relativo de Agua (CRA)}
  \item{y}{Inverso del Potencial Hídrico}
}

\value{Character}
\references{}
\author{Klaus Wolfgang Wiese Acosta}


\examples{
#Set variables
CRA <- c(0.9386,0.9172,0.8791,0.8752,0.8546,0.8411,0.8266,0.8093,0.7759,0.7397,0.7007,0.6828,0.6617)
IPH <- c(2.857,2.500,2.000,1.818,1.666,1.428,1.250,1.086,0.952,0.869,0.781,0.657, 0.606)

#Estimate Parameters from Preassure - Volume Curve
PVCurveMean(1 - CRA, IPH)
}

