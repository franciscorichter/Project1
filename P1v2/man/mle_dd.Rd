\name{mle_dd}
\alias{mle_dd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ MLE for diversity-dependence model
%%  ~~function to do ... ~~
}
\description{ It find the MLE of the dd model given the tree.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mle_dd(n,E,t)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{vector of number of species
%%     ~~Describe \code{x} here~~
}
 \item{E}{binary vector speciation-extinction}
 \item{t}{waiting times}

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
library(DDD)

s3 <- dd_sim(c(0.8,0.1,40),15)
s1 <-L2p(L=s3$L)
p2 <- mle_dd(n=s1$n,E=s1$E,t=s1$t)
c(p2$lambda,p2$beta,p2$mu)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the

