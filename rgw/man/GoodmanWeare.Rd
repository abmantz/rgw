\name{GoodmanWeare}
\alias{GoodmanWeare}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Goodman-Weare Affine-Invariant Sampling
}
\description{
Produces a Monte-Carlo Markov ensemble using the affine-invariant method of Goodman & Weare.
}
\usage{
GoodmanWeare(ensemble, lnpost, Nsteps, current.lnP=NULL,
 mc.cores=getOption("mc.cores", 1L), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ensemble}{
an Nparam*Nwalkers array holding the initial state of the sampler. Nparam is the dimensionality of the parameter space and Nwalkers is the number of positions in the parameter space comprising the ensemble. Nwalkers must be even, and in practice should be *at minimum* twice Nparam.
}
  \item{lnpost}{
function taking a vector of parameter values as input, and returning the log-posterior density.
}
  \item{Nsteps}{
number of iterations to run the sampler.
}
  \item{current.lnP}{
vector holding the log-posterior value corresponding to the initial position of each walker. If not provided, this will be calculated internally.
}
  \item{mc.cores}{
number of cores to use for parallelism.
}
  \item{...}{
additional arguments to pass to lnpost.
}
}
%\details{
%}
\value{
A list containing $ensemble: an array of the same dimensionality as ensemble, containing the position of the walkers after Nsteps iterations of the sampler; and $current.lnP: the  log-posterior density for each walker.
}
\references{
Goodman, J. & Weare, J. (2010, Comm. App. Math. Comp. Sci., 5:6) <DOI:10.2140/camcos.2010.5.65>. This implementation is based on the description given by Foreman-Mackey et al. (2012, arXiv:1202.3665) <DOI:10.1086/670067>.
}
\author{
Adam Mantz
}
\note{
By default, the code will attempt to run in parallel (see the `parallel' package). To prevent this, pass mc.cores=1.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%}
\examples{
# In this example, we'll sample from a simple 2D Gaussian

# Define the log-posterior function
lnP = function(x) sum( dnorm(x, c(0,1), c(pi, exp(0.5)), log=TRUE) )

# Initialize an ensemble of 100 walkers
nwalk = 100
ensemble = array(dim=c(2, nwalk))
ensemble[1,] = rnorm(nwalk, 0, 0.1)
ensemble[2,] = rnorm(nwalk, 1, 0.1)

# Run for a bit
ens2 = GoodmanWeare(ensemble, lnP, 100, mc.cores=1)

# Plot the resulting ensemble
plot(t(ens2$ensemble))
# Compare to a direct draw from the posterior distribution
points(rnorm(nwalk, 0, pi), rnorm(nwalk, 1, exp(0.5)), col=2, pch=3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
