\name{CCMnet}
\alias{CCMnet}
\docType{package}
\title{
Methods to Simulate Networks from the Congruence Class Model.
}
\description{
\pkg{CCMnet} is a collection of methods to simulate networks from the Congruence Class Model (CCM) for networks. For details regarding CCM see the paper by Goyal et al. in the reference section below.

For a list of functions type: help(package='CCMnet')

To cite this package see the citation in \code{citation(package="CCMnet")}. 
}
\details{
The \code{\link{CCMnet_constr}} function simulates networks based on the CCM for networks; see the examples section of the \code{\link{CCMnet_constr}} function for a few demonstrations. 

This package relies on the \pkg{network} package, which allows networks to be represented in \R. Markov Chain Monte
Carlo methods are from the \pkg{ergm} package; \pkg{CCMnet} is partially based on C code from the \pkg{ergm} package (see reference by Handcock et al. below).  The \pkg{SNA} package allows for additional methods relating to social network analysis.  For detailed information regarding these packages go to the \pkg{statnet} website: \url{http://statnet.org}.  This package incorporates an ergm-term coded by Nicole Bohme Carnegie (carnegie@hsph.harvard.edu).
}
\author{Ravi Goyal, with contributions from Mark S. Handcock (ergm R package), David R. Hunter (ergm R package), Carter T. Butts (ergm R package), Steven M. Goodreau (ergm R package), Pavel N. Krivitsky (ergm R package), Martina Morris (ergm R package), Nicole Bohme Carnegie (ergm-user term degmix)
}
\references{
Goyal R, Blitzstein J, and De Gruttola V.  Sampling Networks from Their Posterior Predictive Distribution.  Network Science.  In press.

Handcock M, Hunter D, Butts C, Goodreau S, Krivitsky P, and Morris M (2012). ergm: Fit, Simulate and Diagnose Exponential-Family Models for Networks.  Seattle, WA. Version 3.0-1.  Project home page at http://www.statnet.org.
}