\name{NS_Multinomial}
\alias{NS_Multinomial}
\title{
Calculate network statistic and covariance matrix.
}
\description{
Calculate network statistic and covariance matrix, which is based on a multinomial distribution.  Each unit (either node or edge) in the network is assumed to be sampled from a multinomial distribution based on probabilities associated with the network statistic.
}
\usage{
NS_Multinomial(g,
		   Network_stats, 
		   mean_inflate = 0, 
		   var_inflate = 1,
		   covPattern = NULL
		 )
}
\arguments{
 \item{g}{a network object.}
\item{Network_stats}{Either 'DegreeDist' or 'DegMixing'.}
\item{mean_inflate}{Add small amount to remove zero values from degree mixing matrix entries.}
\item{var_inflate}{Multiply the variance by a constant.  Used to avoid signular covariance matrices.}
\item{covPattern}{Currently not used.}
}
\value{
A list of length 2 containing:
\item{Network Statistic}{Network statistic of the inputted network.}
\item{Covariance}{Covariance matrix for the network statistic; assumes each unit (either node or edge) is sampled from a multinomial distribution based on probabilities derived from the network statistic.}
}
\examples{
g = as.network(rgraph(n=500, m=1, tprob=.01, 
		           mode='graph', diag=FALSE,
                     replace=FALSE, tielist=NULL, 
                     return.as.edgelist=FALSE),
               directed = FALSE)

Prob_Distr_Params=list(NS_Multinomial(g,
				  Network_stats = 'DegreeDist',
                       mean_inflate = .05, 
                       var_inflate = 1.05))

Prob_Distr_Params=list(NS_Multinomial(g,
                       Network_stats = 'DegMixing',
                       mean_inflate = .05, 
                       var_inflate = 1.05))
}