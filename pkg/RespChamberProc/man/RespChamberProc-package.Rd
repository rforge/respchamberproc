\name{RespChamberProc-package}
\alias{RespChamberProc-package}
\alias{RespChamberProc}
\title{Processing data from respiration chambers .}
\description{
Processing data from respiration chambers.
}
\details{

Main method is \code{\link{calcClosedChamberFlux}}, 
to calculate flux and its uncertainty from a time series of measurements.


A first set of functions correct measured CO2 concentrations and fluxes for environmental and chamber properties.
\itemize{
\item{ Dilution with water vapor: \code{\link{corrConcDilution}}  }
\item{ Changing air density with pressure and temperature: \code{\link{corrFluxDensity}}  }
\item{ Leakage of the chamber: \code{\link{corrFluxLeakage}}  }
}%itemize


A second set of functions estimates the flux and its uncertainty from the concentrations. 
\itemize{
\item{ Aggregating estimates: \code{\link{calcClosedChamberFlux}}  }
\item{ Select the subset of data after lagTime: \code{\link{selectDataAfterLag}}  }
\item{ Estimate the initial flux by polynomial regression: \code{\link{regressFluxExp}}  }
\item{ Estimate uncertainty of regression due to leverage of starting and end points: \code{\link{sigmaBootLeverage}}  }
\item{ Convert from mumol CO2/s to g C/day: \code{\link{convert_mumolPers_to_gPerday}}  }
}%itemize

Further misc functionality: \itemize{
\item{ Reading files generated by data-loggers into R: \code{\link{readDat}} and \code{\link{read81x}} }
\item{ Plotting the time series and the resulting best fit: \code{\link{plotResp}} }
\item{ Plotting measurement time versus resulting standard error of the flux estimate: \code{\link{plotDurationUncertainty}} }
}%itemize

Working with a time series of several measurement cycles (see vignette) \itemize{
\item{ identifying contiguous measurement cycles: \code{\link{subsetContiguous}}}
\item{ apply the fitting to all measurement cycles: \code{\link{calcClosedChamberFluxForChunks}}}
\item{ Plotting all time series and fits: \code{\link{plotCampaignConcSeries}} }
}%itemize

}%details

\author{Thomas Wutzler, Oscar Perez Priego}

\keyword{ package }


