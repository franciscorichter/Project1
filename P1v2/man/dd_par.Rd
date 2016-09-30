\name{dd_par}
\alias{dd_par}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Diversity-dependence MLE estimation (parallelized version)
%%  ~~function to do ... ~~
}
\description{ This is the paralelized version of the simulation-estimation framework
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{

u01 <- dd_par(tt=15,mu=0.1,lambda=0.8,K=40,it=10)


}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tt}{Crown time}
  \item{it}{number of iterations}
  \item{mu,lambda,K,etc}{you know this stuff}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
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
library(parallel)
library(foreach)
library(doParallel)

ptm <- proc.time()
u01 <- dd_par(tt=15,mu=0.1,it=10)
proc.time() - ptm
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

