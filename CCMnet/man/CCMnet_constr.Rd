\name{CCMnet_constr}
\alias{CCMnet_constr}
\title{
Congruence Class Model for Networks
}
\description{
Simulate networks from the Congruence Class Model (CCM) for networks.
}
\usage{
CCMnet_constr(Network_stats,
              Prob_Distr, 
              Prob_Distr_Params, 
              samplesize = 5000, 
              burnin=1000, 
              interval=1000,
              statsonly=TRUE,  
              P=NULL,
              population, 
              covPattern = NULL,
              remove_var_last_entry = FALSE,
		   print_info = 0,
              BayesInference = 0, 
              TranNet = NULL, 
              Ia = NULL, 
              Il = NULL, 
              R_times = NULL, 
              beta_a = NULL, 
              beta_l = NULL
		   )
}
\arguments{
 \item{Network_stats}{A vector of network statistics used to generate networks.  Currently implemented network statistics are the following: ``Density'' (Network Density), ``DegreeDist'' (Degree Distribution), ``DegMixing'' (Degree Mixing Matrix), and ``Triangles'' (Number of 3-cycles).}
\item{Prob_Distr}{A vector of probability distributions; one distribution for each network statistic in Network_stats.  Currently implemented probability distributions are: ``Normal'' and ``DirMult''.  DirMult (Dirichlet-multinomial) is currently only implemented for degree distribution.}
\item{Prob_Distr_Params}{A list of lists.  Each list in the list corresponds to a network statistic in Network_stats.  Each list gives the parameters for the associated probability distribution in Prob_Distr. The normal distribution requires a list of length 2 containing a mean vector and covariance matrix.  The Dirichlet-multinomial requires a list of length 1 containing the alpha vector.}
\item{samplesize}{Number of networks simulated.}
\item{burnin}{MCMC burnin.}
\item{interval}{Thinning of MCMC.}
\item{statsonly}{If TRUE, returns only the network statistics of the simulated network, plus the last simulated network.  If FALSE returns the network statistics and simulated networks.}
\item{P}{If P is specified, it will be used as the initial starting network for the MCMC; otherwise, a random graph will be used as the initial starting network.}
\item{population}{Number of nodes in simulated networks.}
\item{covPattern}{Covariate Pattern of each node in network.  Currently only a maximum of two covariate patterns, labeled 1 and 2, are allowed.}
\item{remove_var_last_entry}{Removes the last row and column of the covariance matrix.}
\item{print_info}{Not implemented}
\item{BayesInference}{Not implemented}
\item{TranNet}{Not implemented}
\item{Ia}{Not implemented}
\item{Il}{Not implemented}
\item{R_times}{Not implemented}
\item{beta_a}{Not implemented}
\item{beta_l}{Not implemented}
}
\details{Simulates networks of fixed network size given a probability distribution on network statistics. Uses an algorithm from Goyal et al. and partially based on C code from Handcock et al. (2012).  Also incorporates an ergm-term coded by Nicole Bohme Carnegie (carnegie@hsph.harvard.edu).

Not all combinations of network statistics and probability distributions are currently implemented; additional combinations will be implemented in later releases.  Below we list all implemented combinations.

Density: (Normal)

DegreeDist: (Normal) and (DirMult)

DegMixing: (Normal)

DegMixing and Clustering: (Normal, Normal)

}
\value{
A list of length 2 containing:
\item{Network Statistics}{a matrix of network statistics}
\item{Network(s)}{ List of returned network(s) of type network}
}
\references{
Goyal R, Blitzstein J, and De Gruttola V.  Sampling Networks from Their Posterior Predictive Distribution.  Network Science.  In press.

Handcock M, Hunter D, Butts C, Goodreau S, Krivitsky P, and Morris M (2012). ergm: Fit, Simulate and Diagnose Exponential-Family Models for Networks.  Seattle, WA. Version 3.0-1.  Project home page at http://www.statnet.org.
}
\examples{
set.seed(1)
#
#Example: Density (Normal)
#
CCMnet_Result = CCMnet_constr(Network_stats= c('Density'), 
                              Prob_Distr=c('Normal'), 
                              Prob_Distr_Params=list(list(.04, .0001)), 
                              samplesize = 5000, 
                              burnin=1000, 
                              interval=100,
                              statsonly=TRUE, 
                              P=NULL,
                              population=100, 
                              covPattern = NULL, 
                              remove_var_last_entry = FALSE) 
statsmatrix = CCMnet_Result[[1]]
G_list = CCMnet_Result[[2]]

#Mean Network Density (Simulated Networks)
mean(statsmatrix) 
#Variance of Network Density (Simulated Networks)
var(statsmatrix)

\dontrun{
#
#Example: Degree Distribution (Dirichlet-multinomial)
#
CCMnet_Result = CCMnet_constr(Network_stats='DegreeDist',
                              Prob_Distr='DirMult',
                              Prob_Distr_Params=list(list(c(2,21,15,12))), 
                              samplesize = 10000,
                              burnin=100000, 
                              interval=1000,
                              statsonly=TRUE, 
                              P=NULL,
                              population=500, 
                              covPattern = NULL,
                              remove_var_last_entry = FALSE) 
statsmatrix = CCMnet_Result[[1]]
G_list = CCMnet_Result[[2]]

#Mean Degree Distribution (Simulated Networks)
apply(statsmatrix, 2, mean)
#Variance of Degree Distribution (Simulated Networks)
apply(statsmatrix, 2, var)  

#
#Example: Degree Distribution (Normal)
#
Prob_Distr_Params=list(NS_Multinomial(G_list[[1]],
                                     Network_stats = 'DegreeDist',
                                     mean_inflate = .05, 
                                     var_inflate = 1.05))

CCMnet_Result = CCMnet_constr(Network_stats='DegreeDist',
                              Prob_Distr='Normal',
                              Prob_Distr_Params=Prob_Distr_Params, 
                              samplesize = 50000,
                              burnin=100000, 
                              interval=1000,
                              statsonly=TRUE, 
                              P=NULL,
                              population=500, 
                              covPattern = NULL,
                              remove_var_last_entry = FALSE) 
statsmatrix = CCMnet_Result[[1]]
G_list = CCMnet_Result[[2]]

#Mean Degree Distribution (Simulated Networks)
apply(statsmatrix, 2, mean)
#Variance of Degree Distribution (Simulated Networks)
apply(statsmatrix, 2, var)

#
#Example: Degree Mixing (Normal)
#
Prob_Distr_Params=list(NS_Multinomial(G_list[[1]],
                                     Network_stats = 'DegMixing',
                                     mean_inflate = .05, 
                                     var_inflate = 1.05))

CCMnet_Result = CCMnet_constr(Network_stats='DegMixing',
                              Prob_Distr='Normal',
                              Prob_Distr_Params=Prob_Distr_Params, 
                              samplesize = 50000,
                              burnin=100000, 
                              interval=1000,
                              statsonly=TRUE, 
                              P=NULL,
                              population=500, 
                              covPattern = NULL,
                              remove_var_last_entry = FALSE) 
statsmatrix = CCMnet_Result[[1]]
G = CCMnet_Result[[2]]

#Mean Degree Mixing (Simulated Networks)
apply(statsmatrix, 2, mean)
#Variance of Degree Mixing (Simulated Networks)
apply(statsmatrix, 2, var) 

                   
#
#Example: Degree Mixing and Triangles (Normal, Normal)
#
Prob_Distr_Params=list(NS_Multinomial(G_list[[1]],
                                     Network_stats = 'DegMixing',
                                     mean_inflate = .05, 
                                     var_inflate = 1.05), 
                       list(6,2)) 

CCMnet_Result = CCMnet_constr(Network_stats=c('DegMixing', 'Triangles'),
                              Prob_Distr=c('Normal', 'Normal'),
                              Prob_Distr_Params=Prob_Distr_Params, 
                              samplesize = 50000,
                              burnin=100000, 
                              interval=1000,
                              statsonly=TRUE, 
                              P=NULL,
                              population=500, 
                              covPattern = NULL,
                              remove_var_last_entry = FALSE) 
statsmatrix = CCMnet_Result[[1]]
G = CCMnet_Result[[2]]

#Mean Degree Mixing and Number of Triangles (Simulated Networks)
apply(statsmatrix, 2, mean)
#Variance of Degree Mixing and Number of Triangles (Simulated Networks)
apply(statsmatrix, 2, var)  
}
} 