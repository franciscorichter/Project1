\name{sim_phyl}
\alias{sim_phyl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Phylogenetic tree simulation
%%  ~~function to do ... ~~
}
\description{It simulates a phylogenetic tree given the model parameters
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sim_phyl(ct=15, lambda0=0.8, mu0=0.1, K=40, draw=TRUE, model="dd",printEv=FALSE,seed=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ct}{ Crown time}
  \item{lambda0,mu0}{Speciation and extinction rates parameter (they depends on the model)}
  \item{K}{Carring capacity (For diversity-dependance model only)}
  \item{Draw}{If Draw=TRUE the simulated phylogenetic tree will be ploted}
  \item{Model}{Diversification model. The available models are 'dd' for diversity dependance and 'cr' for constant rate model}
  \item{PrintEv} Prints the events while occurring in the tree
  \item{seed} Seed for the simulation.
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{t} Waiting times of the tree
  \item{E} Binomial vector with 0 for extinctions and 1 for speciations at each waiting time.
  \item{n} Vector with number of species at each waiting time
  \item{newick} Tree in Newick format.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{F. Richter Mendoza
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
st <- sim_phyl2()
plot(st$newick)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

