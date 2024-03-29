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
 save.every=NA, save.file=NA, show.every=NA,
 show.params=1:dim(post)[1], show.walkers=min(dim(post)[2],8),
 show.pch1=1, show.pch2='.', show.pch.switch=500,
 return.lnpost=FALSE, ...)
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
save the accumulated Markov ensemble to disk every time this many
iterations are completed, in a variable called `post`. See `Value`, below.
}
  \item{save.file}{
filename for saving progress.
}
  \item{show.every}{
plot parameter traces so far to the active graphics device periodically.
}
  \item{show.params}{
(sub)set of parameter traces to plot (default is to show all).
}
  \item{show.walkers}{
which walkers to plot traces of (default is first 8).
}
  \item{show.pch1}{
plot symbol to use for short chains.
}
  \item{show.pch2}{
plot symbol to use for long chains.
}
  \item{show.pch.switch}{
chain length that distinguishes "short" and "long" chains for plotting purposes.
}
  \item{return.lnpost}{
whether to return log-posterior values for each sample; see Value.
}
  \item{...}{
additional named arguments to pass to GoodmanWeare or lnpost.
}
}
%\details{
%}
\value{
If return.lnpost==FALSE, an array of the same dimensionality as post, storing the position of the walkers in post[,,i] every thin iterations. Otherwise, a list containing that array as $post, as well as an Nwalkers*Nsteps array storing the corresponding log-posterior values as $lnP. The log-posterior values $lnP[,1], corresponding with the starting ensemble positions $post[,,1], will always be NA.
}
\references{
See also help for rgw::GoodmanWeare.
}
\author{
Adam Mantz
}
\note{
By default, the code will attempt to run in parallel (see the `parallel'
package). To prevent this, pass mc.cores=1.

If traces are being plotted (show.every not NA),
par(mfrow=c(length(show.params), 1)) is called on the current graphics
device. If the device is not large enough to show all of the traces,
this will cause a crash.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%}
\examples{
# In this example, we'll sample from a simple 2D Gaussian.
# (This is the same example as used in GoodmanWeare.)

# Define the log-posterior function
lnP = function(x) sum( dnorm(x, c(0,1), c(pi, exp(0.5)), log=TRUE) )

# Initialize an ensemble of 100 walkers. We'll take 100 steps, saving the
# ensemble after each.
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
