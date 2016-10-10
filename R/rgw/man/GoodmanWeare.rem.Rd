\name{GoodmanWeare.rem}
\alias{GoodmanWeare.rem}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Goodman-Weare Affine-Invariant Sampling
}
\description{
Produces a Monte-Carlo Markov ensemble using the affine-invariant method of Goodman & Weare, saving progress periodically.
}
\usage{
GoodmanWeare.rem(post, lnpost, thin=1, mention.every=NA,
 save.every=NA, save.file=NA, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{post}{
an Nparam*Nwalkers*Nsteps array. post[,,1] should hold the initial state of the sampler (see help for GoodmanWeare). Checkpoints and the return value will have the same shape, with subsequent layers post[,,i] holding the ensemble state at later iterations.
}
  \item{lnpost}{
function taking a vector of parameter values as input, and returning the log-posterior density.
}
  \item{thin}{
thinning factor for saving the results.
}
  \item{mention.every}{
print a message to the console every time this many iterations are completed.
}
  \item{save.every}{
save the accumulated Markov ensemble to disk every time this many iterations are completed.
}
  \item{save.file}{
filename for saving progress.
}
  \item{...}{
additional arguments to pass to GoodmanWeare or lnpost.
}
}
%\details{
%}
\value{
An array of the same dimensionality as post, storing the position of the walkers in post[,,i] every thin iterations.
}
\references{
See help for GoodmanWeare.
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
# In this example, we'll sample from a simple 2D Gaussian.
# (This is the same example as used in GoodmanWeare.)

# Define the log-posterior function
lnP = function(x) sum( dnorm(x, c(0,1), c(pi, exp(0.5)), log=TRUE) )

# Initialize an ensemble of 100 walkers. We'll take 100 steps, saving the ensemble after each.
nwalk = 100
post = array(NA, dim=c(2, nwalk, 101))
post[1,,1] = rnorm(nwalk, 0, 0.1)
post[2,,1] = rnorm(nwalk, 1, 0.1)

# Run
post = GoodmanWeare.rem(post, lnP, mc.cores=1)

# Plot the final ensemble
plot(post[1,,101], post[2,,101])
# Look at the trace of each parameter for one of the walkers.
plot(post[1,1,])
plot(post[2,1,])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
